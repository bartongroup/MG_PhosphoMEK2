
make_metadata <- function(samples, conditions) {
  samples %>% 
    separate(sample, c("condition", "replicate"), sep = "-", remove = FALSE) %>% 
    mutate(
      condition = factor(condition, levels = conditions),
      replicate = factor(replicate)
    ) %>% 
    arrange(condition, replicate)
}

read_mq <- function(file, data_cols, measure_cols, id_cols, filt, meta) {
  mr <- meta %>% 
    full_join(measure_cols, by = "reporter")
  n2n <- set_names(data_cols$raw_name, data_cols$name)
  raw <- read_tsv(file, col_select = c(data_cols$raw_name, mr$column_name), show_col_types = FALSE) %>% 
    rename(all_of(n2n)) %>% 
    mutate(id = as.character(id))
  dat <- raw %>% 
    filter(rlang::eval_tidy(rlang::parse_expr(filt))) %>% 
    select(all_of(id_cols), all_of(measure_cols$column_name)) %>% 
    pivot_longer(-all_of(id_cols), names_to = "column_name", values_to = "value") %>% 
    mutate(value = na_if(value, 0)) %>% 
    left_join(mr, by = "column_name") %>% 
    select(id, multi, sample, value) %>% 
    mutate(multi = multi %>% as.character())
  info <- raw %>% select(-all_of(mr$column_name))
  
  set <- list(
    info = info,
    dat = dat,
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
    left_join(rep_cols, by = "column_name")
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

# At least one non-missing value in all replicates in at least one condition.
get_expressed_ids <- function(set) {
  set$dat %>% 
    left_join(set$metadata, by = "sample") %>% 
    group_by(id, multi, condition) %>% 
    summarise(n_tot = n(), n_good = length(na.omit(value))) %>% 
    ungroup() %>% 
    group_by(id, multi) %>% 
    summarise(n_good_conditions = sum(n_tot == n_good)) %>% 
    filter(n_good_conditions > 0) %>% 
    select(id, multi)
}

set_comparison <- function(pho, pep, pro) {
  good_pho <- get_expressed_ids(pho)
  good_pep <- get_expressed_ids(pep)
  good_pro <- get_expressed_ids(pro)
  
  pho2pep <- pho$info %>% 
    select(phospho_id = id, peptide_id = peptide_ids) %>% 
    separate_rows(peptide_id, sep = ";") %>% 
    filter(phospho_id %in% good_pho$id & peptide_id %in% good_pep$id) %>% 
    distinct()
  
  pho2pro <- pho$info %>% 
    select(phospho_id = id, protein_id = protein_ids) %>% 
    separate_rows(protein_id, sep = ";") %>% 
    filter(phospho_id %in% good_pho$id & protein_id %in% good_pro$id) %>% 
    distinct()
  
  tb <- good_pho %>%
    select(phospho_id = id, multi) %>% 
    full_join(pho2pep, by = "phospho_id") %>%
    full_join(pho2pro, by = "phospho_id")
  
  list(
    `Phospho site` = unique(tb$phospho_id),
    `With peptide` = unique(tb %>% filter(!is.na(peptide_id)) %>% pull(phospho_id)),
    `With protein` = unique(tb %>% filter(!is.na(protein_id)) %>% pull(phospho_id))
  )
}

# match all phosphosites vs all proteins, mark expressed as "good"
pho_pro_match <- function(pho, pro) {
  good_pho <- get_expressed_ids(pho) %>% 
    add_column(good_pho = 1)
  good_pro <- get_expressed_ids(pro) %>% 
    select(-multi) %>% 
    add_column(good_pro = 1)
  
  pho2pro <- pho$info %>% 
    select(phospho_id = id, protein_id = protein_ids) %>% 
    separate_rows(protein_id, sep = ";") %>% 
    distinct()
  
  pro2pho <- pro$info %>% 
    select(protein_id = id, phospho_id = phospho_ids) %>% 
    separate_rows(phospho_id, sep = ";") %>% 
    distinct()
  
  bind_rows(pho2pro, pro2pho) %>%
    distinct() %>% 
    left_join(good_pho, by = c("phospho_id" = "id")) %>%
    left_join(good_pro, by = c("protein_id" = "id")) %>%
    filter(!(is.na(good_pho) & is.na(good_pro))) %>% 
    mutate(
      good_pho = !is.na(good_pho),
      good_pro = !is.na(good_pro)
    )
}


get_phospho_genes <- function(pho) {
  pho$info %>% 
    select(gene_name) %>% 
    drop_na() %>% 
    separate_rows(gene_name, sep = ";") %>% 
    distinct() %>% 
    pull(gene_name)
}

# number of quantified features
# quantified in all replicates in at least one condition
n_quantified <- function(set) {
  get_expressed_ids(set) %>% 
    length()
}


# quantified phospho sites, peptides and proteins
q_numbers <- function(pho, pep, pro) {
  map_dfr(list(pho, pep, pro), function(x) {tibble(n = n_quantified(x))}) %>% 
    add_column(set = c("Phospho", "Peptide", "Protein"), .before = 1)
}


#' RAS procedure for CONSTANd normalization (not exported)
#'
#' @param K Input matrix
#' @param max.iter Maximum number of iterations
#' @param eps Convergence limit
#'
#' @return Matrix normalized so row and column means equal 1/n (n - number of
#'   columns)
RAS <- function(K, max.iter = 50, eps = 1e-5) {
  n <- ncol(K)
  m <- nrow(K)
  row_names <- rownames(K)
  col_names <- colnames(K)
  
  # ignore rows with only NAs
  good.rows <- which(rowSums(!is.na(K)) > 0)
  K <- K[good.rows, ]
  
  cnt <- 1
  repeat {
    row.mult <- 1 / (n * rowMeans(K, na.rm = TRUE))
    K <- K * row.mult
    err1 <- 0.5 * sum(abs(colMeans(K, na.rm = TRUE) - 1/n))
    col.mult <- 1 / (n * colMeans(K, na.rm = TRUE))
    K <- t(t(K) * col.mult)
    err2 <- 0.5 * sum(abs(rowMeans(K, na.rm = TRUE) - 1/n))
    cnt <- cnt + 1
    if (cnt > max.iter || (err1 < eps && err2 < eps)) break
  }
  
  # reconstruct full table
  KF <- matrix(NA, nrow = m, ncol = n)
  KF[good.rows, ] <- K
  
  colnames(KF) <- col_names
  rownames(KF) <- row_names
  
  return(KF)
}


normalise_constand <- function(set) {
  tab <- dat2mat(set$dat, "value")
  tab_norm <- RAS(tab)
  dn <- tab_norm %>% 
    as_tibble(rownames = "mid") %>% 
    separate(mid, c("id", "multi"), sep = "-") %>% 
    pivot_longer(-c(id, multi), names_to = "sample", values_to = "value_constand")
  set$dat <- set$dat %>% 
    left_join(dn, by = c("id", "multi", "sample"))
  set
}

normalise_to_median <- function(set) {
  med <- set$dat %>% 
    group_by(sample) %>% 
    summarise(M = median(value, na.rm = TRUE)) %>% 
    mutate(M = M / mean(M))
  set$dat <- set$dat %>% 
    left_join(med, by = "sample") %>% 
    mutate(value_med = value / M) %>% 
    select(-M)
  set
}


normalise_to_proteins <- function(pho, pro) {
  # here we ignore a handful of phospho sites that a linked to multiple protein groups
  pho$phospho2prot <- pho$info %>% 
    select(id, protein_id = protein_ids) %>% 
    filter(!str_detect(protein_id, ";"))
  # mean protein abundance across conditions
  mp <- pro$dat %>% 
    drop_na() %>% 
    left_join(pro$metadata, by = "sample") %>% 
    group_by(id, condition) %>% 
    summarise(prot_mean = mean(value)) %>% 
    rename(protein_id = id)
  pho$dat <- pho$dat %>% 
    left_join(pho$phospho2prot, by = "id") %>% 
    left_join(select(pho$metadata, sample, condition), by = "sample") %>% 
    left_join(mp, by = c("protein_id", "condition")) %>% 
    left_join(select(pro$dat, protein_id = id, sample, prot = value), by = c("protein_id", "sample")) %>% 
    mutate(
      value_prot = value / prot,
      value_prot_mean = value / prot_mean
    ) %>% 
    select(-c(protein_id, condition, prot, prot_mean))
  pho
}


dat2mat <- function(d, what) {
  d %>% 
    select(id, multi, sample, val = !!what) %>% 
    unite(mid, c(id, multi), sep = "-") %>% 
    drop_na() %>% 
    pivot_wider(id_cols = mid, names_from = sample, values_from = val) %>% 
    column_to_rownames("mid") %>% 
    as.matrix()
}

# Protein per phosphosite statistics. How many proteins a phosphosite is matched
# to? 
protein_count <- function(tab_pho_pro) {
  tab_pho_pro %>%
    filter(good_pho) %>%  
    group_by(good_pro, phospho_id, .drop = FALSE) %>%
    summarise(n_prot = n()) %>%
    ungroup() %>% 
    group_by(good_pro, n_prot) %>%
    summarise(count = n()) %>% 
    arrange(good_pro, n_prot) %>% 
    ungroup() %>% 
    mutate(`Protein detected` = if_else(good_pro, "Yes", "No"), .before = good_pro) %>% 
    select(-good_pro) %>% 
    rename(`N matched proteins` = n_prot)
}


detect_duplicates <- function(set) {
  dups <- set$dat %>% 
    select(id, multi, sample, value) %>%
    pivot_wider(id_cols = c(id, multi), values_from = value, names_from = sample) %>%
    drop_na() %>% 
    unite("values", all_of(set$metadata$sample)) %>% 
    group_by(values) %>%
    mutate(group_id = cur_group_id(), n = n()) %>%
    ungroup() %>%
    filter(n > 1) %>% 
    arrange(group_id)
}

convert_sequences <- function(win, wcov, pos) {
  sqw <- win %>% str_split("") %>% unlist()
  pw <- wcov %>% str_split("") %>% unlist()
  sqw[pw == "P"] %>%
    str_c(collapse = "") %>% 
    mark_position(pos)
}


get_phospho_info <- function(pep, pho, pho_del) {
  pho_info <- pho$info %>% 
    right_join(pho_del, by = "id") %>% 
    select(phospho_id = id, multi, peptide_ids, protein, gene_name, localization_prob, amino_acid, position)
  
  peptide_ids <- pho_info$peptide_ids %>% 
    str_split(pattern = ";") %>% 
    unlist()
  
  pep_info <- pep$info %>% 
    filter(id %in% peptide_ids) %>% 
    select(peptide_id = id, phospho_ids, sequence, start_position, end_position) %>% 
    # if one phospho site is in multiple peptides, select peptide with longest sequence
    mutate(length = nchar(sequence)) %>% 
    arrange(desc(length), sequence) %>% 
    slice(1)
  
  bind_cols(pho_info, pep_info) %>% 
    mutate(position_in_peptide = position - start_position + 1) %>% 
    mutate(sequence_pho = mark_position(sequence, position_in_peptide))
}


duplicate_example <- function(pep, pho, dup, gr) {
  pho_sel <- dup %>% 
    filter(group_id == gr) %>% 
    select(id, multi)
  get_phospho_info(pep, pho, pho_sel) %>% 
    select(phospho_id, multi, protein, gene_name, localization_prob, amino_acid, position_in_peptide, position, sequence)


}