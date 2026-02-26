targets_main <- function() {

  biomart <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    tar_target(bm_genes, biomart_fetch_genes(mart)),
    #tar_target(go_terms,  bm_fetch_go(mart, phospho_genes)),
    #tar_target(gs_terms,  bm_fetch_go(mart, phospho_genes, slim = TRUE)),
    #tar_target(re_terms, fetch_reactome(mart, phospho_genes)),
    #tar_target(kg_terms, get_kegg(species = KEGG_SPECIES, bm_genes = bm_genes)),
    #tar_target(all_terms, list(go = go_terms, gs = gs_terms, re = re_terms, kg = kg_terms)),
    tar_target(terms, get_functional_terms_sym()),
    tar_target(fterms, prepare_terms_fenr(terms, phospho_genes))
  )

  read_data <- list(
    tar_target(metadata, make_metadata(SAMPLE_REPORTER, CONDITIONS$CONDITION)),
    tar_target(proteins, read_mq(PROTEINS_FILE, PROTEINS_DATA_COLUMNS, PROTEINS_MEASURE_COLUMNS, PROTEINS_ID_COLUMNS, PROTEINS_FILTER, metadata)),
    tar_target(peptides, read_mq(PEPTIDES_FILE, PEPTIDES_DATA_COLUMNS, PEPTIDES_MEASURE_COLUMNS, PEPTIDES_ID_COLUMNS, PEPTIDES_FILTER, metadata)),
    tar_target(phospho,
      read_mq(PHOSPHO_FILE, PHOSPHO_DATA_COLUMNS, PHOSPHO_MEASURE_COLUMNS, PHOSPHO_ID_COLUMNS, PHOSPHO_FILTER, metadata) %>%
      normalise_to_proteins(proteins) %>%
      deduplicate()
    ),
    tar_target(phospho_rep, read_phospho_reporters(PHOSPHO_FILE))
  )

  stats <- list(
    tar_target(phospho_rep_names, phospho_rep %>% pull(column_name) %>% unique()),
    tar_target(quants, q_numbers(phospho, peptides, proteins)),
    tar_target(upset_pho_pep_pro, set_comparison(phospho, peptides, proteins)),
    tar_target(n_good_phospho, prepare_phospho_counts(phospho, loc_prob_limit = 0.75) %>% pull(id) %>% unique() %>% length()),
    tar_target(phospho_dup_example, duplicate_example(peptides, phospho, 74)),
    tar_target(n_dup_removed, n_duplicates(phospho$duplicates))
  )

  proteins <- list(
    tar_target(peptide_de, limma_de(peptides, info_cols = KEEP_PEPTIDES_COLUMNS, what = "value_med", log_scale = TRUE)),
    tar_target(protein_de, limma_de(proteins, info_cols = KEEP_PROTEINS_COLUMNS, what = "value_med", log_scale = TRUE)),
    tar_target(fig_volcano_prot, plot_volcano(protein_de) + scale_y_continuous(expand = c(0, 0), limits = c(0, 10))),

    tar_target(tab_pho_pro, pho_pro_match(phospho, proteins)),
    tar_target(tabs_prot_count, protein_count(tab_pho_pro))
  )

  map_normalisations <- tar_map(
    values = NORMALISATIONS,
    names = NAME,

    tar_target(phospho_de, limma_de(phospho, info_cols = KEEP_PHOSPHO_COLUMNS, what = WHAT, log_scale = LOG, loc_prob_limit = 0.75)),
    tar_target(fig_volcano, plot_volcano(phospho_de, logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT) + scale_y_continuous(expand = c(0, 0), limits = c(0, 10))),
    tar_target(fig_ma, plot_ma(phospho_de, logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT)),
    tar_target(fig_up_down, plot_up_down(phospho_de, logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT)),
    tar_target(fig_sample_dist, plot_sample_distirbutions(phospho, WHAT, log_scale = LOG, ncol = 5)),
    tar_target(fig_pho_distmat, plot_distance_matrix(phospho, WHAT)),
    tar_target(fig_pho_clustering, plot_clustering(phospho, WHAT)),

    tar_target(de_sites, phospho_de %>% filter(FDR < FDR_LIMIT & abs(logFC) >= LOGFC_LIMIT) %>% select(id, multi)),
    tar_target(fig_pho_vs_pro, plot_pho_vs_prot(phospho, proteins, de_sites)),
    tar_target(tab_de, make_de_table(phospho_de, phospho$info))
  )

  selections <- list(
    tar_target(phospho_genes, get_phospho_genes(phospho)),
    tar_target(genes_de, make_de_genes(phospho_de_median, fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT))
  )

  stringdb <- tar_map(
    values = tibble::tibble(direction = c("up", "down")),
    names = direction,
    tar_target(stringr_map, run_stringdb(genes_de[[direction]]$protein, species = TAXONOMY_ID)),
    tar_target(png_strigr, plot_stringdb_clusters(stringr_map, paste0("fig/sdb_", direction, ".png")))
  )

  figures <- list(
    tar_target(upset_reporters, upset_phospho_orders_overlap(phospho_rep)),
    tar_target(fig_phorep_1, plot_phospho_orders(phospho_rep, 1)),
    tar_target(fig_prot_norm_problem, plot_phospho_norm(phospho, proteins, pho_id = 23999, pho_multi = 1)),
    tar_target(fig_pho_per_pep, plot_phospho_per_peptide(peptides, phospho)),
    tar_target(fig_pho_pro_cor, plot_pho_prot_cor(phospho, proteins)),
    tar_target(fig_pho_dups, plot_duplications(phospho$duplicates)),

    tar_target(fig_de_heatmap, plot_de_heatmap(phospho, de_sites_median, what = "value_med"))
  )

  figures_pairs <- tar_map(
    values = NORM_COND,
    names = NAME,
    tar_target(fig_pairs, plot_replicate_pairs(phospho, CONDITION, WHAT, log_scale = LOG))
  )

  shiny <- list(
    tar_target(sav_shiny_de, shiny_data_de(phospho, peptides, proteins, phospho_de_median, bm_genes, fterms) %>% write_rds("shiny/data_de.rds", compress = "xz"))
  )

  tables <- list(
    tar_target(sav_de, tab_de_median %>% write_tsv("tab/phospho_de.tsv"))
  )

  pd_comparison <- list(
    tar_target(pd, read_and_process_pd_data(PROTEOME_DISCOVERER_FILE, metadata)),
    tar_target(pd_mq, map_mq_pd(phospho$info, pd$info)),
    tar_target(pd_mq_comp, plot_pd_mq_comparison(pd, phospho, peptides, pd_mq, c("Q9UPT8:S1269", "Q9UPT8:S1275")))
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
    tables,
    pd_comparison
  )

}
