
make_metadata <- function(samples, conditions) {
  samples %>%
    separate(sample, c("condition", "replicate"), sep = "-", remove = FALSE) %>%
    mutate(
      condition = factor(condition, levels = conditions),
      replicate = factor(replicate)
    ) %>%
    arrange(condition, replicate)
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



compare_pd_mq <- function(pd, mq, pd_mq) {
  pd_mq <- drop_na(pd_mq)
  abu <- pd$abu %>%
    left_join(pd$info, by = "peptide_id") %>%
    select(pd_id = peptide_id, pd_multi = multi, sample, abundance, abundance_norm)
  dat <- mq$dat %>%
    select(mq_id = id, mq_multi = multi, sample, value, value_med)
  
  d <- dat %>% 
    drop_na() %>% 
    inner_join(pd_mq, by = c("mq_id")) %>% 
    inner_join(abu, by = c("sample", "pd_id"))
}


plot_pd_mq_comparison <- function(pd, mq, mq_pep, pd_mq, prot_pos) {
  pm <- pd_mq %>% 
    filter(prot_position %in% prot_pos)
  mq_ids <- pm$mq_id %>% unique()
  pd_ids <- pm$pd_id %>% unique()
  
  mq_info <- mq$info %>% 
    filter(id %in% mq_ids)
  pd_info <- pd$info %>% 
    filter(peptide_id %in% pd_ids)
  
  pep_ids <- mq_info %>% 
    select(peptide_ids) %>% 
    separate_rows(peptide_ids, sep = ";") %>% 
    pull(peptide_ids) %>% 
    as.integer() %>% 
    unique()
  pep_info <- mq_pep$info %>% 
    filter(id %in% pep_ids)
  
  d <- compare_pd_mq(pd, mq, pm)
  
  g0 <- ggplot(d) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0)
    ) +
    scale_colour_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    labs(x = NULL)
  g1 <- g0 +
    geom_point(aes(x = sample, y = value, colour = as.factor(mq_id), shape = as.factor(mq_multi))) +
    scale_shape_manual(values = c(15, 16)) +
    labs(colour = "Site ID", shape = "Multiplicity", title = "MaxQuant")
  g2 <- g0 +
    geom_point(aes(x = sample, y = abundance, colour = as.factor(pd_id))) +
    labs(colour = "Peptide ID", title = "Proteome Discoverer")
  g <- plot_grid(g1, g2, nrow = 1, align = "h")
  
  list(
    mq_info = mq_info,
    pd_info = pd_info,
    pep_info = pep_info,
    plot = g
  )
}