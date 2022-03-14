targets_main <- function() {
  
  biomart <- list(
    tar_target(mart, useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version="105")),
    tar_target(bm_genes, biomart_fetch_genes(mart)),
    tar_target(bm_go, bm_fetch_go(mart, detected_genes)),
    tar_target(bm_go_slim, bm_fetch_go(mart, detected_genes, slim=TRUE)),
    tar_target(reactome, fetch_reactome(mart, detected_genes)),
    tar_target(kegg, get_kegg(species = "mmu", bm_genes = bm_genes))
  )
  
  read_data <- list(
    tar_target(metadata, make_metadata(SAMPLE_REPORTER, CONDITIONS$condition)),
    tar_target(proteins, read_mq(PROTEINS_FILE, PROTEINS_DATA_COLUMNS, PROTEINS_MEASURE_COLUMNS, PROTEINS_ID_COLUMNS, metadata)),
    #tar_target(peptides, read_mq(PEPTIDES_FILE, PEPTIDES_DATA_COLUMNS, PEPTIDES_MEASURE_COLUMNS, PEPTIDES_ID_COLUMNS, metadata)),
    tar_target(phospho, read_mq(PHOSPHO_FILE, PHOSPHO_DATA_COLUMNS, PHOSPHO_MEASURE_COLUMNS, PHOSPHO_ID_COLUMNS, metadata) %>% 
      normalise_to_proteins(proteins)),
    tar_target(phospho_rep, read_phospho_reporters(PHOSPHO_FILE))
  )
  
  stats <- list(
    tar_target(phospho_rep_names, phospho_rep %>% pull(column_name) %>% unique()),
    #tar_target(quants, q_numbers(phospho, peptides, proteins)),
    tar_target(detected_genes, get_detected_genes(phospho)),
    #tar_target(upset_pho_pep_pro, set_comparison(phospho, peptides, proteins)),
    tar_target(n_good_phospho, prepare_phospho_counts(phospho, tab = "tab_med", loc_prob_limit = 0.95) %>% nrow())
  )
  
  de <- list(
    tar_target(phospho_de, limma_de(phospho, info_cols=KEEP_PHOSPHO_COLUMNS, tab = "tab_med", log_scale = TRUE, loc_prob_limit=0.95)),
    tar_target(phospho_n_de, limma_de(phospho, info_cols=KEEP_PHOSPHO_COLUMNS, tab = "tab_prot", log_scale = TRUE, loc_prob_limit=0.95)),
    #tar_target(peptide_de, limma_de(peptides, info_cols=KEEP_PEPTIDES_COLUMNS, tab = "tab_med", log_scale = TRUE)),
    tar_target(protein_de, limma_de(proteins, info_cols=KEEP_PROTEINS_COLUMNS, tab = "tab_med", log_scale = TRUE)),
    
    tar_target(fig_volcano, plot_volcano(phospho_de, logfc.limit = 1)),
    tar_target(fig_ma, plot_ma(phospho_de, logfc.limit = 1)),
    tar_target(fig_up_down, plot_up_down(phospho_de, fc.limit=1)),
    tar_target(upset_fc1, de_list(phospho_de, "contrast", logfc.limit = 1))
  )
  
  fgsea <- list(
    tar_target(fg_go, fgsea_cache(phospho_de, bm_go, "cache/fgsea_go.rds")),
    tar_target(fg_gs, fgsea_cache(phospho_de, bm_go_slim, "cache/fgsea_gs.rds")),
    tar_target(fg_re, fgsea_cache(phospho_de, reactome, "cache/fgsea_re.rds")),
    tar_target(fg_kg, fgsea_cache(phospho_de, kegg, "cache/fgsea_kg.rds"))
  )
  
 
  figures <- list(
    tar_target(upset_reporters, upset_phospho_orders_overlap(phospho_rep)),
    tar_target(fig_phorep_1, plot_phospho_orders(phospho_rep, 1)),
    tar_target(fig_sample_dist_med, plot_sample_distirbutions(phospho, "value_med", log_scale = TRUE)),
    tar_target(fig_sample_dist_constand, plot_sample_distirbutions(phospho, "value_constand", log_scale = FALSE))
  )
  
  figures_pairs <- tar_map(
    values = CONDITIONS,
    tar_target(fig_pair_med, plot_replicate_pairs(phospho, condition, what = "value_med", log_scale = TRUE)),
    tar_target(fig_pair_constand, plot_replicate_pairs(phospho, condition, what = "value_constand", log_scale = FALSE))
  )
  
  shiny <- list(
    tar_target(shiny_de, shiny_data_de(phospho, peptides, proteins, phospho_de, bm_genes, bm_go, reactome, kegg)),
    tar_target(sav_shiny_de, write_rds(shiny_de, "shiny/data_de.rds", compress = "xz"))
  )
  
  tables <- list(
    tar_target(sav_de, phospho_de %>% mutate_if(is.numeric, ~signif(.x, 4)) %>% write_tsv("tab/phospho_de.tsv"))
  )
  
  c(
    biomart,
    read_data,
    de,
    #fgsea,
    stats,
    figures,
    figures_pairs
    #shiny,
    #tables
  )
  
}