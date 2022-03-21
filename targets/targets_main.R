targets_main <- function() {
  
  biomart <- list(
    tar_target(mart, useEnsembl(biomart="ensembl", dataset=ENSEMBL_DATASET, version=ENSEMBL_VERSION)),
    tar_target(bm_genes, biomart_fetch_genes(mart)),
    tar_target(all_terms, list(
      go = bm_fetch_go(mart, phospho_genes),
      gs = bm_fetch_go(mart, phospho_genes, slim=TRUE),
      re = fetch_reactome(mart, phospho_genes),
      kg = get_kegg(species = KEGG_SPECIES, bm_genes = bm_genes))
    )
  )
  
  read_data <- list(
    tar_target(metadata, make_metadata(SAMPLE_REPORTER, CONDITIONS$CONDITION)),
    tar_target(proteins, read_mq(PROTEINS_FILE, PROTEINS_DATA_COLUMNS, PROTEINS_MEASURE_COLUMNS, PROTEINS_ID_COLUMNS, PROTEINS_FILTER, metadata)),
    tar_target(peptides, read_mq(PEPTIDES_FILE, PEPTIDES_DATA_COLUMNS, PEPTIDES_MEASURE_COLUMNS, PEPTIDES_ID_COLUMNS, PEPTIDES_FILTER,metadata)),
    tar_target(phospho, read_mq(PHOSPHO_FILE, PHOSPHO_DATA_COLUMNS, PHOSPHO_MEASURE_COLUMNS, PHOSPHO_ID_COLUMNS, PHOSPHO_FILTER, metadata) %>% 
      normalise_to_proteins(proteins)),
    tar_target(phospho_rep, read_phospho_reporters(PHOSPHO_FILE))
  )
  
  stats <- list(
    tar_target(phospho_rep_names, phospho_rep %>% pull(column_name) %>% unique()),
    tar_target(quants, q_numbers(phospho, peptides, proteins)),
    tar_target(upset_pho_pep_pro, set_comparison(phospho, peptides, proteins)),
    tar_target(n_good_phospho, prepare_phospho_counts(phospho, loc_prob_limit = 0.95) %>% pull(id) %>% unique() %>% length())
  )
  
  proteins <- list(
    tar_target(peptide_de, limma_de(peptides, info_cols=KEEP_PEPTIDES_COLUMNS, what = "value_med", log_scale = TRUE)),
    tar_target(protein_de, limma_de(proteins, info_cols=KEEP_PROTEINS_COLUMNS, what = "value_med", log_scale = TRUE)),
    tar_target(fig_volcano_prot, plot_volcano(protein_de)),
    
    tar_target(tab_pho_pro, pho_pro_match(phospho, proteins)),
    tar_target(tabs_prot_count, protein_count(tab_pho_pro))
  )
  
  map_normalisations <- tar_map(
    values = NORMALISATIONS,
    names = NAME,
  
    tar_target(phospho_de, limma_de(phospho, info_cols=KEEP_PHOSPHO_COLUMNS, what = WHAT, log_scale = LOG, loc_prob_limit = 0.75)),
    tar_target(fig_volcano, plot_volcano(phospho_de, logfc.limit = LOGFC_LIMIT, fdr.limit = FDR_LIMIT)),
    tar_target(fig_ma, plot_ma(phospho_de, logfc.limit = LOGFC_LIMIT, fdr.limit = FDR_LIMIT)),
    tar_target(fig_up_down, plot_up_down(phospho_de, logfc.limit = LOGFC_LIMIT, fdr.limit = FDR_LIMIT)),
    tar_target(fig_sample_dist, plot_sample_distirbutions(phospho, WHAT, log_scale = LOG, ncol=5)),
    tar_target(fig_pho_distmat, plot_distance_matrix(phospho, WHAT)),
    tar_target(fig_pho_clustering, plot_clustering(phospho, WHAT)),
    
    tar_target(de_sites, phospho_de %>% filter(FDR < FDR_LIMIT & abs(logFC) >= LOGFC_LIMIT) %>% select(id, multi)),
    tar_target(fig_pho_vs_pro, plot_pho_vs_prot(phospho, proteins, de_sites)),
    tar_target(tab_de, make_de_table(phospho_de, phospho$info))
  )
  
  selections <- list(
    tar_target(phospho_genes, get_phospho_genes(phospho)),
    tar_target(genes_de, make_de_genes(phospho_de_median, fdr.limit=FDR_LIMIT, logfc.limit=LOGFC_LIMIT))
  )
  
  stringdb <- tar_map(
    values = tibble::tibble(direction = c("up", "down")),
    names = direction,
    tar_target(stringr_map, run_stringdb(genes_de[[direction]]$protein, species=TAXONOMY_ID)),
    tar_target(png_strigr, plot_stringdb_clusters(stringr_map, paste0("fig/sdb_", direction, ".png")))
  )
  
  figures <- list(
    tar_target(upset_reporters, upset_phospho_orders_overlap(phospho_rep)),
    tar_target(fig_phorep_1, plot_phospho_orders(phospho_rep, 1)),
    tar_target(fig_prot_norm_problem, plot_phospho_norm(phospho, proteins, pho_id="23999", pho_multi="1")),
    tar_target(fig_pho_per_pep, plot_phospho_per_peptide(peptides, phospho)),
    tar_target(fig_pho_pro_cor, plot_pho_prot_cor(phospho, proteins)),
    
    tar_target(fig_de_heatmap, plot_de_heatmap(phospho, de_sites_median, what = "value_med"))
  )
  
  figures_pairs <- tar_map(
    values = NORM_COND,
    names = NAME,
    tar_target(fig_pairs, plot_replicate_pairs(phospho, CONDITION, WHAT, log_scale = LOG))
  )

  shiny <- list(
    tar_target(sav_shiny_de, shiny_data_de(phospho, peptides, proteins, phospho_de_median, bm_genes, all_terms) %>% write_rds("shiny/data_de.rds", compress = "xz"))
  )
  
  tables <- list(
    tar_target(sav_de, tab_de_median %>% write_tsv("tab/phospho_de.tsv"))
  )
  
  c(
    biomart,
    read_data,
    proteins,
    map_normalisations,
    stats,
    selections,
    stringdb,
    figures,
    figures_pairs,
    shiny,
    tables
  )
  
}