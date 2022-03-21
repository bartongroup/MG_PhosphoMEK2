okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


plot_sample_distirbutions <- function(set, what="value_norm", bins=50, log_scale = FALSE, ncol=3) {
  d <- set$dat %>% 
    drop_na() %>% 
    mutate(val = get(what))
  if(log_scale) d$val = log10(d$val)
  dm <- d %>% 
    group_by(sample) %>% 
    summarise(M = median(val))
  g <- ggplot(d, aes(x=val, y=..density..)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(bins = bins, fill = "grey50") +
    #geom_vline(data = dm, aes(xintercept = M)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = what, y = NULL) +
    facet_wrap(~sample, ncol=ncol)
  plot_grid(g)
}

plot_replicate_pairs <- function(set, cond, what="value_med", log_scale = FALSE) {
  dat <- set$dat %>% 
    mutate(val = get(what)) %>% 
    select(id, multi, sample, val) %>% 
    left_join(set$metadata, by="sample") %>% 
    filter(condition == cond) %>% 
    drop_na() %>% 
    select(id, multi, sample, val)
  if(log_scale) dat$val <- log10(dat$val)
  samples <- dat$sample %>% unique()
  g <- expand_grid(x=samples, y=samples) %>% 
    filter(x != y) %>%
    left_join(dat, by=c("x" = "sample")) %>%
    left_join(dat, by=c("id", "multi", "y" = "sample")) %>% 
  ggplot(aes(x=val.x, y=val.y)) +
    theme_bw() + 
    theme(panel.grid = element_blank(), legend.position = "none") +
    #geom_point(size = 0.3, colour="grey30", alpha=0.3) +
    geom_hex(bins=200) +
    geom_abline(slope = 1, intercept = 0, colour="brown") +
    facet_grid(y ~ x) +
    labs(x = NULL, y = NULL)
  plot_grid(g)
}

plot_distance_matrix <- function(set, what = "value_med", text.size=10) {
  tab <- dat2mat(set$dat, what) %>% 
    log10()

  cor(tab, use="complete.obs") %>% 
    as_tibble(rownames = "sample") %>%
    pivot_longer(-sample) %>% 
    mutate(sample = factor(sample, levels=set$metadata$sample)) %>% 
    mutate(name = factor(name, levels=set$metadata$sample)) %>% 
    ggplot(aes(x=sample, y=name)) +
    geom_tile(aes(fill=value)) +
    scale_fill_viridis(option="cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=text.size),
      axis.text.y = element_text(size=text.size)
    ) +
    labs(x=NULL, y=NULL, fill="Correlation")
}


plot_clustering <- function(set, what = "value_med", text.size=10) {
  tab <- dat2mat(set$dat, what) %>% 
    log10()
  
  corr.mat <- cor(tab, use="complete.obs")
  dis <- as.dist(1 - corr.mat)  # dissimilarity matrix
  hc <- hclust(dis)
  dendr <- ggdendro::dendro_data(hc)
  seg <- ggdendro::segment(dendr)
  meta <- set$metadata %>% mutate(sample = as.character(sample))
  labs <- left_join(dendr$labels %>% mutate(label = as.character(label)), meta, by=c("label"="sample")) %>% 
    mutate(colour = okabe_ito_palette[as_factor(condition)])
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size=text.size, colour=labs$colour),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(size=0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data=seg, aes_(x=~x, y=~y, xend=~xend, yend=~yend)) +
    scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values=okabe_ito_palette) +
    labs(x=NULL, y="Distance")
}

plot_volcano <- function(res, fc="logFC", p="PValue", fdr="FDR", groupvar="contrast",
                         fdr.limit=0.05, logfc.limit=0, point.size=0.5, point.alpha=0.5) {
  res %>%
    mutate(
      x = get(fc),
      y = get(p),
      sig = get(fdr) < fdr.limit & abs(get(fc)) > logfc.limit
    ) %>%
    ggplot(aes(x=x, y=-log10(y), colour=sig)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_vline(xintercept = 0, colour = "brown") +
    geom_point(size=point.size, alpha=point.alpha) +
    scale_colour_manual(values=c("grey70", "black")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_grid(as.formula(glue::glue(". ~ {groupvar}"))) +
    theme(legend.position = "none") +
    labs(x=fc, y=glue::glue("-log10({p})"))
}


plot_ma <- function(res, fc="logFC", sm="AveExpr", fdr = "FDR", groupvar = "contrast",
                    fdr.limit=0.05, logfc.limit=0, with.labels=FALSE, label.size=3, id="gene_name",
                    point.size=0.5, point.alpha=0.5, sel_genes=NULL) {
  res <- res %>%
    mutate(sig = get(fdr) < fdr.limit & abs(get(fc)) > logfc.limit)
  g <- ggplot(res, aes_string(x=sm, y=fc, colour="sig")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 0, colour = "brown") +
    geom_point(size=point.size, alpha=point.alpha) +
    scale_colour_manual(values=c("grey70", "black")) +
    facet_grid(as.formula(glue::glue(". ~ {groupvar}"))) +
    theme(legend.position = "none")
  if(!is.null(sel_genes)) {
    g <- g + geom_point(data=filter(res, res$gene_name %in% sel_genes), colour="red", size=point.size*1.5)
  }
  if(with.labels) {
    g <- g + geom_text_repel(data=filter(res, sig), aes_string(x=sm, y=fc, label=id), size=label.size)
  }
  g
}


plot_up_down <- function(res, fdr.limit=0.05, logfc.limit=1, groupvar = "contrast") {
  res %>%
    filter(FDR < fdr.limit & abs(logFC) > logfc.limit) %>%
    mutate(direction = if_else(logFC > 0, "up", "down")) %>%
    mutate(gr = !!sym(groupvar)) %>% 
    group_by(gr, direction) %>%
    tally() %>%
    mutate(n = if_else(direction == "down", -n, n)) %>%
    ggplot(aes(x = gr, y = n, fill=direction)) +
    theme_bw() +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette) +
    labs(x=NULL, y="Num significant genes")
}


mark_position <- function(s, pos) {
  s <- tolower(s) 
  str_sub(s, start=pos, end=pos) <- toupper(str_sub(s, start=pos, end=pos))
  s
}

plot_phospho_intensities <- function(pho, phospho_id, what="value_med", log_scale=TRUE, tit=NULL) {
  pho_dat <- pho$dat %>% 
    filter(id == phospho_id) %>% 
    mutate(val = get(what)) %>% 
    select(id, sample, val) %>% 
    drop_na() %>% 
    left_join(pho$metadata, by="sample")
  if(nrow(pho_dat) == 0) return(NULL)
  
  if(log_scale) pho_dat$val <- log10(pho_dat$val)
  
  pho_info <- pho$info %>% 
    filter(id == phospho_id)
  
  # extract peptide sequence
  if(is.null(tit)) {
    sqw <- pho_info$sequence_window %>% str_split("") %>% unlist()
    pw <- pho_info$peptide_window_coverage %>% str_split("") %>% unlist()
    sq <- sqw[pw == "P"] %>%
      str_c(collapse = "") %>% 
      mark_position(pho_info$position_in_peptide)
    tit <- glue::glue("{sq}:{pho_info$position} {pho_info$gene_name}")
  }
  
  ggplot(pho_dat, aes(x=condition, y=val, colour=replicate)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    scale_colour_manual(values = okabe_ito_palette) +
    geom_point() +
    labs(x=NULL, y=what, title=tit)
}


plot_peptide <- function(pep, pho, peptide_id) {
  pep_info <- pep$info %>% 
    filter(id == peptide_id)
  phospho_ids <- get_phospho_ids(pep, peptide_id)
  
  map(phospho_ids, function(ph_id) {
    plot_phospho(pho, ph_id)
  }) %>% 
    compact() %>% 
    plot_grid(plotlist = ., nrow = 1)
}


plot_peptide_seq <- function(pep, pho, de, peptide_ids) {
  pep_info <- pep$info %>% 
    filter(id %in% peptide_ids) %>% 
    select(peptide_id = id, phospho_id = phospho_ids, sequence, start_position, end_position) %>% 
    separate_rows(phospho_id, sep=";") %>% 
    # if one phospho site is in multiple peptides, select peptide with longest sequence
    mutate(length = nchar(sequence)) %>% 
    arrange(desc(length)) %>% 
    group_by(phospho_id) %>% 
    slice(1) %>% 
    ungroup()

  pho_info <- pho$info %>% 
    select(phospho_id = id, protein, gene_name, localization_prob, amino_acid, position)

  pho_de <- de %>% 
    select(phospho_id = id, logFC, AveExpr, FDR, contrast)
  
  d <- pep_info %>% 
    left_join(pho_info, by = "phospho_id") %>% 
    left_join(pho_de, by = "phospho_id") %>% 
    mutate(position_in_peptide = position - start_position + 1) %>% 
    drop_na()
  mx <- max(abs(na.omit(d$logFC)))
  
  d %>% 
    group_split(peptide_id) %>% 
    map(function(w) {
      dl <- w$position_in_peptide
      sq <- first(w$sequence) %>% str_split("") %>% unlist()
      title <- paste0(first(w$protein), ":", first(w$gene_name), " ", first(w$start_position), "-", first(w$end_position))
      dp <- tibble(
        aa = sq,
        pos = seq_along(aa)
      ) %>% 
        left_join(w, by = c("pos" = "position_in_peptide"))
      n <- nrow(dp)

    ggplot(dp, aes(x=as_factor(pos), y=logFC)) +
      theme_void() +
      theme(
        axis.ticks.y = element_line(size=0.5),
        axis.ticks.length = unit(.1, "cm"),
        axis.text.y = element_text(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(fill = alpha("lightgoldenrod1", 0.3), colour=NA)
      ) +
      geom_col(position = position_dodge(), aes(fill = contrast, colour = factor(FDR < 0.05, levels=c(FALSE, TRUE)))) +
      geom_hline(yintercept = 0, size=2, colour="grey50", alpha=0.2) +
      geom_text(aes(x=pos, y=0, label=aa), size=5) +
      scale_x_discrete(breaks = dp$pos, labels = dp$aa, drop=FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(-mx, mx), breaks = scales::breaks_width(2)) +
      scale_fill_manual(values = okabe_ito_palette, drop=FALSE) +
      scale_colour_manual(values = c("grey70", "black"), drop=FALSE) +
      ggtitle(title)
    }) %>% 
    plot_grid(plotlist = ., scale=0.9)
}



plot_phospho_per_peptide <- function(pep, pho) {
  # only quantified phosphosites (across all samples)
  pho_quant <- pho$dat %>% 
    group_by(id) %>% 
    summarise(n_tot = n(), n_good = (length(na.omit(value)))) %>% 
    filter(n_good == n_tot) %>% 
    pull(id)
  
  g <- pep$info %>% 
    select(id, phospho_ids) %>% 
    separate_rows(phospho_ids, sep=";") %>% 
    filter(phospho_ids %in% pho_quant) %>% 
    group_by(id) %>% 
    tally() %>% 
    group_by(cnt = n) %>% 
    tally() %>% 
  ggplot(aes(x=cnt, y=n)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_x_continuous(breaks = 1:30) +
    labs(x="Number of phosphosites", y="Peptide count")
  plot_grid(g)
}


plot_full_protein <- function(de, pro, pro_id, cntr) {
  this_pro <- pro$info %>% 
    filter(id == pro_id) 
  this_pro %>% 
    separate_rows(phospho_ids, sep = ";") %>% 
    select(id = phospho_ids, sequence_length) %>% 
    left_join(de) %>% 
    filter(!is.na(logFC) & contrast == cntr) %>% 
    ggplot(aes(x=position, y=logFC, colour=FDR<0.01)) + 
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    geom_point() +
    geom_segment(aes(xend=position, yend=0)) +
    geom_hline(yintercept = 0) +
    scale_colour_manual(values = c("grey80", "black")) +
    scale_x_continuous(expand = c(0,0), limits = c(0, this_pro$sequence_length)) +
    labs(x = "Position", y=expression(log[2]~FC), title=glue::glue("{this_pro$protein} : {this_pro$gene_name}"))
}



upset_phospho_orders_overlap <- function(pr) {
  orders <- unique(pr$order)
  d <- pr %>%
    mutate(value = na_if(value, 0)) %>%
    pivot_wider(id_cols = c(id, reporter), names_from = order, values_from = value) %>% 
    unite("idr", c(id, reporter))
  map(orders, function(ord) {
    d %>%
      select(idr, all_of(ord)) %>% 
      drop_na() %>% 
      pull(idr)
  }) %>% 
    set_names(orders)
}


plot_phospho_orders <- function(pr, rep=1) {
  pr %>%
    filter(reporter == rep) %>% 
    mutate(value = na_if(value, 0)) %>%
    mutate(value = log10(value)) %>% 
    pivot_wider(id_cols = c(id, reporter), names_from = order, values_from = value) %>% 
    select(-c(id, reporter)) %>% 
  ggpairs()
}


# A look at normalisation to protein
plot_phospho_norm <- function(pho, pro, pho_id, pho_multi) {
  g <- pho$dat %>% 
    filter(id == pho_id & multi == pho_multi) %>% 
    left_join(pho$phospho2prot, by = "id") %>% 
    left_join(pro$dat %>% select(protein_id = id, sample, prot = value), by = c("protein_id", "sample")) %>% 
    select(-c(id, multi, value_constand, protein_id)) %>% 
    rename(
      `Protein raw` = prot,
      `Phospho raw` = value,
      `Phospho to median` = value_med,
      `Phospho to protein` = value_prot,
      `Phospho to mean protein` = value_prot_mean,
    ) %>% 
    pivot_longer(-sample) %>% 
    left_join(pho$metadata, by = "sample") %>% 
  ggplot(aes(x = condition, y = value, colour = condition)) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    #geom_segment(aes(xend = condition, yend = 0), colour = "grey") +
    geom_point(size=3) +
    scale_colour_manual(values = okabe_ito_palette) +
    #scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    facet_wrap(~name, scales = "free_y", nrow=1) +
    labs(x=NULL, y="Intensity or ratio")
  plot_grid(g)
}


plot_pho_vs_prot <- function(pho, pro, sites, ncol=6) {
  pho$dat %>%
    right_join(sites, by = c("id", "multi")) %>% 
    left_join(pho$metadata, by = "sample") %>% 
    left_join(pho$phospho2prot, by = "id") %>% 
    left_join(pro$dat %>% select(protein_id = id, sample, prot = value_med), by = c("protein_id", "sample")) %>% 
    mutate(prot = replace_na(prot, 1)) %>% 
    unite(mid, c("id", "multi"), sep = "-") %>% 
  ggplot(aes(x = log10(value_med), y = log10(prot), fill = condition)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point(shape = 21, size = 2, colour = "grey20") +
    scale_fill_manual(values = c("#FFFFFF", "#000000")) +
    labs(x = expression(log[10]~phosphosite), y = expression(log[10]~protein)) +
    facet_wrap(~mid, scales = "free", ncol=ncol)
}

plot_pho_prot_cor <- function(pho, pro, nrow=2) {
  d <- pho$dat %>%
    drop_na() %>% 
    left_join(pho$phospho2prot, by = "id") %>% 
    full_join(select(pro$dat, id, sample, prot = value), by = c("sample", "protein_id" = "id")) %>% 
    mutate(pho = log10(value), prot = log10(prot)) %>% 
    select(id, sample, pho, prot)
  dc <- d %>% 
    group_by(sample) %>% 
    summarise(cor = cor(pho, prot))
  g <- d %>% 
    ggplot(aes(x = prot, y = pho)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point(size = 0.05) +
    facet_wrap(~ sample, nrow=nrow) +
    labs(x = expression(log[10]~protein), y = expression(log[10]~phosphosite))
  plot_grid(g)
}


plot_de_heatmap <- function(set, sites, what = "value_med", max.scale=NULL) {
  if(nrow(sites) == 0) return(NULL)
  X <- set$dat %>%
    mutate(val = get(what)) %>% 
    right_join(sites, by = c("id", "multi")) %>% 
    mutate(sample = factor(sample, levels=set$metadata$sample)) %>% 
    group_by(id, multi) %>%
    mutate(M = mean(val), logfc = log2(val / M)) %>%
    pivot_wider(id_cols = c(id, multi), names_from = sample, values_from = logfc) %>%
    left_join(set$info %>% select(id, gene_name), by="id") %>% 
    unite("gid", gene_name, id, multi, sep="-") %>% 
    column_to_rownames("gid")
  ggheatmap(X, legend.name=expression(log[2]~FC), with.y.text = TRUE)
}

