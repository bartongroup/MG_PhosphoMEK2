
make_metadata <- function(samples, conditions) {
  samples %>% 
    separate(sample, c("condition", "replicate"), sep = "-", remove = FALSE) %>% 
    mutate(
      condition = factor(condition, levels = conditions),
      replicate = factor(replicate)
    ) %>% 
    arrange(condition, replicate)
}

read_mq <- function(file, data_cols, measure_cols, id_cols, meta) {
  mr <- meta %>% 
    left_join(measure_cols, by="reporter")
  raw <- read_tsv(file, col_select = c(data_cols$raw_name, mr$column_name), show_col_types = FALSE) %>% 
    set_names(c(data_cols$name, mr$sample)) %>% 
    mutate(id = as.character(id)) %>% 
    filter(!str_detect(protein, "^(CON_|REV_)")) # remove reverse and contaminants
  dat <- raw %>% 
    select(all_of(id_cols), all_of(meta$sample)) %>%
    pivot_longer(-all_of(id_cols), names_to = "sample", values_to = "value") %>% 
    mutate(sample = factor(sample, levels = meta$sample)) %>% 
    mutate(value = na_if(value, 0))
  tab <- raw[, meta$sample] %>% as.matrix()
  rownames(tab) <- raw$id
  tab[tab == 0] <- NA
  info <- raw %>% select(-all_of(meta$sample))
  
  set <- list(
    info = info,
    dat = dat,
    tab = tab,
    metadata = meta
  )
  
  set <- normalise_to_median(set)
  set <- normalise_constand(set)
  
  set
}

# For testing
read_phospho_reporters <- function(file) {
  rep_cols <- expand_grid(reporter = REPORTERS, order = 1:3) %>% 
    mutate(column_name = glue::glue("Reporter intensity corrected {reporter}___{order}") %>% as.character())
  read_tsv(file, col_select = c("id", rep_cols$column_name), show_col_types = FALSE) %>% 
    pivot_longer(-id, names_to = "column_name") %>% 
    left_join(rep_cols, by="column_name")
}


select_full_rows <- function(tab) {
  tab[rowSums(is.na(tab)) == 0, ]
}

drop_empty_rows <- function(tab) {
  tab[rowSums(!is.na(tab)) > 0, ]
}


get_peptide_ids <- function(pho, phospho_ids) {
  pho$info %>% 
    filter(id %in% phospho_ids) %>% 
    pull(peptide_ids) %>% 
    str_split(";") %>% 
    unlist() %>% 
    unique()
}

get_phospho_ids <- function(pep, peptide_ids) {
  pep$info %>% 
    filter(id %in% peptide_ids) %>% 
    pull(phospho_ids) %>% 
    str_split(";") %>% 
    unlist() %>% 
    unique()
}


get_expressed_ids <- function(set) {
  set$dat %>% 
    left_join(set$metadata, by="sample") %>% 
    group_by(id, condition) %>% 
    summarise(n_good = sum(!is.na(value_med))) %>% 
    ungroup() %>% 
    group_by(id) %>% 
    summarise(good = sum(n_good) > 0) %>% 
    arrange(as.numeric(id)) %>% 
    filter(good) %>% 
    pull(id)
}

set_comparison <- function(pho, pep, pro) {
  good_pho <- get_expressed_ids(pho)
  good_pep <- get_expressed_ids(pep)
  good_pro <- get_expressed_ids(pro)
  
  pho2pep <- pho$info %>% 
    select(phospho_id = id, peptide_id = peptide_ids) %>% 
    separate_rows(peptide_id, sep=";") %>% 
    filter(phospho_id %in% good_pho & peptide_id %in% good_pep) %>% 
    distinct()
  
  pho2pro <- pho$info %>% 
    select(phospho_id = id, protein_id = protein_ids) %>% 
    separate_rows(protein_id, sep=";") %>% 
    filter(phospho_id %in% good_pho & protein_id %in% good_pro) %>% 
    distinct()
  
  tb <- tibble(phospho_id = good_pho) %>%
    full_join(pho2pep) %>%
    full_join(pho2pro)
  
  list(
    `Phospho site` = unique(tb$phospho_id),
    `With peptide` = unique(tb %>% filter(!is.na(peptide_id)) %>% pull(phospho_id)),
    `With protein` = unique(tb %>% filter(!is.na(protein_id)) %>% pull(phospho_id))
  )
}


get_detected_genes <- function(pho) {
  pho$dat %>% 
    select(id, value) %>% 
    drop_na() %>% 
    group_by(id) %>% 
    tally() %>% 
    left_join(pho$info, by="id") %>% 
    select(gene_name) %>% 
    separate_rows(gene_name, sep=";") %>% 
    distinct() %>% 
    pull(gene_name)
}

# number of quantified features
# quantified in all replicates in at least one condition
n_quantified <- function(set) {
  set$dat %>% 
    left_join(set$metadata, by="sample") %>% 
    group_by(id, condition) %>% 
    summarise(n_tot = n(), n_good = length(na.omit(value))) %>% 
    ungroup() %>% 
    group_by(id) %>% 
    summarise(n_good_conditions = sum(n_tot == n_good)) %>% 
    filter(n_good_conditions > 0) %>% 
    nrow()
}


# quantified phospho sites, peptides and proteins
q_numbers <- function(pho, pep, pro) {
  map_dfr(list(pho, pep, pro), function(x) {tibble(n = n_quantified(x))}) %>% 
    add_column(set = c("Phospho", "Peptide", "Protein"), .before=1)
}


#' RAS procedure for CONSTANd normalization (not exported)
#'
#' @param K Input matrix
#' @param max.iter Maximum number of iterations
#' @param eps Convergence limit
#'
#' @return Matrix normalized so row and column means equal 1/n (n - number of
#'   columns)
RAS <- function(K, max.iter=50, eps=1e-5) {
  n <- ncol(K)
  m <- nrow(K)
  row_names <- rownames(K)
  col_names <- colnames(K)
  
  # ignore rows with only NAs
  good.rows <- which(rowSums(!is.na(K)) > 0)
  K <- K[good.rows, ]
  
  cnt <- 1
  repeat {
    row.mult <- 1 / (n * rowMeans(K, na.rm=TRUE))
    K <- K * row.mult
    err1 <- 0.5 * sum(abs(colMeans(K, na.rm=TRUE) - 1/n))
    col.mult <- 1 / (n * colMeans(K, na.rm=TRUE))
    K <- t(t(K) * col.mult)
    err2 <- 0.5 * sum(abs(rowMeans(K, na.rm=TRUE) - 1/n))
    cnt <- cnt + 1
    if(cnt > max.iter || (err1 < eps && err2 < eps)) break
  }
  
  # reconstruct full table
  KF <- matrix(NA, nrow=m, ncol=n)
  KF[good.rows, ] <- K
  
  colnames(KF) <- col_names
  rownames(KF) <- row_names
  
  return(KF)
}


normalise_constand <- function(set) {
  tab_norm <- RAS(set$tab)
  dn <- tab_norm %>% 
    as_tibble(rownames = "id") %>% 
    pivot_longer(-id, names_to = "sample", values_to = "value_constand")
  set$dat <- set$dat %>% 
    left_join(dn, by = c("id", "sample"))
  set$tab_constand <- tab_norm
  set
}

normalise_to_median <- function(set) {
  normf <- apply(set$tab, 2, function(x) median(x, na.rm=TRUE))
  normf <- normf / mean(normf)
  tab_norm <- t(t(set$tab) / normf)
  dn <- tab_norm %>% 
    as_tibble(rownames = "id") %>% 
    pivot_longer(-id, names_to = "sample", values_to = "value_med")
  set$dat <- set$dat %>% 
    left_join(dn, by = c("id", "sample"))
  set$tab_med <- tab_norm
  set
}


normalise_to_proteins <- function(pho, pro) {
  # here we ignore a handful of phospho sites that a linked to multiple protein groups
  p2p <- pho$info %>% 
    select(id, protein_id = protein_ids) %>% 
    filter(!str_detect(protein_id, ";"))
  dat <- pho$dat %>% 
    left_join(p2p, by="id") %>% 
    left_join(select(pro$dat, protein_id = id, sample, prot = value, prot_norm = value_med), by=c("protein_id", "sample")) %>% 
    mutate(ratio_norm = value_med / prot_norm)
  tab_prot <- dat %>% 
    pivot_wider(id_cols = id, names_from = sample, values_from = ratio_norm) %>% 
    column_to_rownames("id") %>% 
    as.matrix()
  list(
    info = pho$info,
    phospho2prot = p2p,
    dat = dat,
    tab = pho$tab,
    tab_med = pho$tab_med,
    tab_constand = pho$tab_constand,
    tab_prot = tab_prot,
    metadata = pho$metadata
  )
}
