okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


sh_plot_volcano <- function(d, alpha = 0.05, title = NULL) {
  sres <- d %>% filter(FDR <= alpha)
  g <- ggplot(d, aes(x, y)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(size = 0.2) +
    geom_vline(xintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P), title = title) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
  if (nrow(sres) > 0) {
    g <- g + geom_hline(yintercept = -log10(max(sres$PValue)), colour = "red", linetype = "dashed", alpha = 0.2)
  }
  g
}

sh_plot_ma <- function(d, alpha = 0.05, title = NULL) {
  ggplot(d, aes(x, y)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(data = d[d$FDR > 0.05,], size = 0.1, colour = "grey50") +
    geom_point(data = d[d$FDR <= 0.05,], size = 0.2, colour = "black") +
    geom_hline(yintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~CPM), y = expression(log[2]~FC), title = title) 
}


sh_functional_enrichment <- function(genes_all, genes_sel, term_data, gene2name = NULL,
                                 min_count = 2, sig_limit = 0.05) {
  
  gene2term <- term_data$gene2term %>% mutate(gene_name = toupper(gene_name))
  term_info <- term_data$terms
  
  # select only terms represented in our gene set
  gene2term <- gene2term %>% filter(gene_name %in% genes_all)
  
  # all terms present in the selection
  terms <- gene2term %>% 
    filter(gene_name %in% genes_sel) %>% 
    pull(term_id) %>% 
    unique()
  
  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)
  
  # empty line for missing terms
  na_term <- term_info %>% slice(1) %>% mutate_all(~NA)
  
  res <- map_dfr(terms, function(term) {
    info <- term_info %>% filter(term_id == term)
    # returns NAs if no term found
    if (nrow(info) == 0) info <- na_term %>% mutate(term_id = term)
    
    # all genes with the term
    tgenes <- gene2term %>% filter(term_id == term) %>% pull(gene_name)
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)
    
    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)
    
    expected <- nuni * Nsel / Nuni
    fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow = 2)
    ft <- fisher.test(fish, alternative = "greater")
    p <- as.numeric(ft$p.value)
    
    if (!is.null(gene2name)) tgenes_sel <- gene2name[tgenes_sel] %>% unname()
    
    bind_cols(
      info,
      tibble(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        ids = paste(tgenes_sel, collapse = ","),
        P = p
      )
    )
  }) %>% 
    mutate(P = p.adjust(P, method = "BH")) %>% 
    filter(sel >= min_count & P <= sig_limit) %>% 
    arrange(desc(enrich)) %>% 
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))
  
  res
}

all_table <- function(d) {
  d %>% 
    select(id, logFC, FDR, protein, gene_name) %>% 
    mutate_at(c("logFC", "FDR"), ~signif (.x, 3)) %>% 
    DT::datatable(class = 'cell-border strip hover', selection = "single", rownames = FALSE) 
}


sh_mark_position <- function(s, pos) {
  s <- tolower(s) 
  str_sub(s, start = pos, end = pos) <- toupper(str_sub(s, start = pos, end = pos))
  s
}


sh_get_phospho_info <- function(pep, pho, phospho_id) {
  pho_info <- pho$info %>% 
    filter(id == phospho_id) %>% 
    select(phospho_id = id, peptide_ids, protein, gene_name, localization_prob, amino_acid, position)
  
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
    mutate(sequence_pho = sh_mark_position(sequence, position_in_peptide))
}


sh_get_all_phosphos <- function(pep, pho, peptide_id) {
  pep_info <- pep$info %>% 
    filter(id == peptide_id) %>% 
    select(peptide_id = id, phospho_ids, sequence, start_position, end_position) %>% 
    separate_rows(phospho_ids, sep = ";") %>% 
    mutate(phospho_ids = as.integer(phospho_ids)) %>% 
    rename(phospho_id = phospho_ids) %>% 
    left_join(pho$info %>% select(phospho_id = id, peptide_ids, protein, gene_name, localization_prob, amino_acid, position), by = "phospho_id") %>% 
    mutate(position_in_peptide = position - start_position + 1) %>% 
    mutate(sequence_pho = sh_mark_position(sequence, position_in_peptide))
}


sh_plot_pepseq <- function(pep, pho, de, pho_id) {
  merged <- sh_get_phospho_info(pep, pho, pho_id)
  pep_pho <- sh_get_all_phosphos(pep, pho, merged$peptide_id)
  
  pho_de <- de %>% 
    select(phospho_id = id, logFC, AveExpr, FDR, contrast)
  
  w <- pep_pho %>% 
    left_join(pho_de, by = "phospho_id") 
  
  mx <- max(abs(na.omit(w$logFC)))
  dl <- w$position_in_peptide
  sq <- first(w$sequence) %>% str_split("") %>% unlist()
  title <- paste0(first(w$protein), ":", first(w$gene_name), " ", first(w$start_position), "-", first(w$end_position))
  dp <- tibble(
    aa = sq,
    pos = seq_along(aa)
  ) %>% 
    left_join(w, by = c("pos" = "position_in_peptide")) %>% 
    mutate(this_one = phospho_id == pho_id)
  txt <- dp %>%
    select(pos, aa, this_one) %>% 
    distinct() %>% 
    mutate(this_one = replace_na(this_one, FALSE))
  n <- nrow(dp)
  
  ggplot(dp, aes(x = as_factor(pos), y = logFC)) +
    theme_void() +
    theme(
      axis.ticks.y = element_line(size = 0.5),
      axis.ticks.length = unit(.1, "cm"),
      axis.text.y = element_text(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.background = element_rect(fill = alpha("lightgoldenrod1", 0.3), colour = NA)
    ) +
    geom_col(position = position_dodge(), aes(fill = contrast, colour = factor(FDR < 0.05, levels = c(FALSE, TRUE)))) +
    geom_hline(yintercept = 0, size = 2, colour = "grey50", alpha = 0.2) +
    geom_text(data = txt %>% filter(!this_one), aes(x = pos, y = 0, label = aa), size = 5, colour = "black") +
    geom_text(data = txt %>% filter(this_one), aes(x = pos, y = 0, label = aa), size = 8, colour = "red") +
    scale_x_discrete(breaks = dp$pos, labels = dp$aa, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(-mx, mx), breaks = scales::breaks_width(2)) +
    scale_fill_manual(values = okabe_ito_palette, drop = FALSE, na.translate = FALSE) +
    scale_colour_manual(values = c("grey70", "black"), drop = FALSE) +
    guides(color = "none", fill = guide_legend(title = NULL))
}


sh_plot_full_protein <- function(de, pro, pro_id, pho_id, cntr) {
  this_pro <- pro$info %>% 
    filter(id %in% pro_id) 
  if (nrow(this_pro) > 1) this_pro = this_pro[1, ]
  d <- this_pro %>% 
    separate_rows(phospho_ids, sep = ";") %>% 
    mutate(phospho_ids = as.integer(phospho_ids)) %>% 
    select(id = phospho_ids, sequence_length) %>% 
    left_join(de, by = "id") %>% 
    filter(!is.na(logFC) & contrast == cntr)
 ggplot(d, aes(x = position, y = logFC, colour = FDR < 0.01)) + 
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    geom_segment(aes(xend = position, yend = 0)) +
    geom_point(size = 2) +
    geom_point(data = d %>% filter(id == pho_id), colour = "red", size = 3) +
    geom_hline(yintercept = 0) +
    scale_colour_manual(values = c("grey80", "black")) +
    scale_x_continuous(expand = c(0,0), limits = c(0, this_pro$sequence_length)) +
    labs(x = "Position", y = expression(log[2]~FC), title = NULL)
}



sh_plot_intensities <- function(set, pid, what = "value_med", log_scale = TRUE, tit = NULL) {
  dat <- set$dat %>% 
    filter(id == pid) %>% 
    mutate(val = get(what)) %>% 
    select(id, sample, val) %>% 
    drop_na() %>% 
    left_join(set$metadata, by = "sample")
  if (nrow(dat) == 0) return(NULL)
  
  if (log_scale) {
    dat$val <- log10(dat$val)
    ylab <- expression(log[10]~Intensity)
  } else {
    dat$val <- dat$val / 1e6
    ylab <- expression(Intensity~x~10^6)
  }
  mx <- max(dat$val) * 1.05
  
  g <- ggplot(dat, aes(x = condition, y = val, colour = replicate, group = replicate)) + 
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    scale_colour_manual(values = okabe_ito_palette) +
    geom_beeswarm(size = 3, cex = 1.5) +
    labs(x = NULL, y = ylab, title = tit)
  if (!log_scale) g <- g + scale_y_continuous(expand = c(0, 0), limits = c(0, mx))
  g
}



