prepare_phospho_counts <- function(pho, loc_prob_limit=0.75) {
  sel <- pho$info %>% 
    filter(localization_prob > loc_prob_limit) %>% 
    pull(id)
  pho$dat %>%
    filter(id %in% sel)
}


limma_de <- function(set, formula="~ 0 + condition", sig.level=0.05, info_cols=NULL, what = "value_med", log_scale=TRUE, loc_prob_limit=0.75) {
  if("localization_prob" %in% colnames(set$info)) {
    d <- prepare_phospho_counts(set, loc_prob_limit)
  } else {
    d <- set$dat
  }
  
  X <- dat2mat(d, what)
  
  if(log_scale) X <- log2(X)
  meta <- set$metadata
  conditions <- levels(meta$condition)
  
  design <- model.matrix(as.formula(formula), meta)
  colnames(design) <- conditions

  ctrs <- expand_grid(x=as_factor(conditions), y=as_factor(conditions)) %>%
    filter(as.integer(x) < as.integer(y)) %>%
    unite(contrast, c(y, x), sep="-") %>%
    pull(contrast)
  contrast_mat <- makeContrasts(contrasts=ctrs, levels=design)
  
  fit <- lmFit(X, design) %>% 
    contrasts.fit(contrasts = contrast_mat) %>% 
    eBayes()
  
  
  res <- map_dfr(ctrs, function(cf) {
    topTable(fit, coef=cf, number=1e6, sort.by="none") %>%
      as_tibble(rownames = "mid") %>%
      separate(mid, c("id", "multi"), sep = "-") %>% 
      mutate(contrast = cf) %>%
      rename(FDR = adj.P.Val, PValue = P.Value) %>%
      select(-c(t, B))
  }) %>% 
    drop_na() %>% 
    mutate(contrast = factor(contrast, levels = ctrs))
  
  if(!is.null(info_cols)) {
    info <- set$info %>% 
      select(c("id", all_of(info_cols)))
    res <- res %>% 
      left_join(info, by="id")
  }
  res
}




de_list <- function(res, group_var, fdr="FDR", logfc="logFC", logfc.limit=0, fdr.limit=0.05, name=NULL, split_up_down=FALSE) {
  d <- res %>% 
    filter(!!sym(fdr) < fdr.limit & abs(!!sym(logfc)) >= logfc.limit) %>% 
    mutate(group = !!sym(group_var))
  if(split_up_down) {
    d <- d %>% 
      mutate(direction = if_else(!!sym(logfc) > 0, "up", "down")) %>% 
      unite(group, c(group, direction), sep=":")
  }  
  d <- d %>% 
    select(group, id) %>%
    group_by(group)
  kname <- ifelse(is.null(name), "", paste0(name, ":"))
  ks <- paste0(kname, group_keys(d)[[1]])
  d %>%
    distinct() %>%
    group_map(~pull(.x, id)) %>%
    set_names(ks)
}


pull_proteins <- function(des) {
  des %>% 
    select(protein, gene_name) %>% 
    distinct() %>% 
    arrange(gene_name) %>% 
    filter(!is.na(gene_name))
}

make_de_genes <- function(de, fdr.limit=0.05, logfc.limit=1) {
  list(
    up = pull_proteins(de %>% filter(FDR < fdr.limit & logFC >= logfc.limit)),
    down = pull_proteins(de %>% filter(FDR < fdr.limit & logFC <= -logfc.limit))
  )
}

make_de_table <- function(de, info, fdr.limit, logfc.limit) {
  de %>%
    #filter(FDR < fdr.limit & abs(logFC) >= logfc.limit) %>% 
    select(id, multi, logFC, PValue, FDR) %>% 
    left_join(info, by = "id") %>% 
    select(id, multi, logFC, FDR, proteins, protein_names, gene_name, aa = amino_acid, in_pep = position_in_peptide, in_prot = position, charge) %>% 
    mutate(across(c(logFC, FDR), ~signif(.x, 3)))
}