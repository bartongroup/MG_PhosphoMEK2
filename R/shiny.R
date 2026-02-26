shiny_data_de <- function(phospho, peptides, proteins, phospho_de, bm_genes, fterms) {
  pho2pep <- phospho$info %>%
    select(id, peptide_ids) %>%
    separate_rows(peptide_ids, sep = ";") %>% 
    mutate(peptide_ids = as.integer(peptide_ids))
  pho2pro <- phospho$info %>%
    select(id, protein_ids) %>%
    separate_rows(protein_ids, sep = ";") %>% 
    mutate(protein_ids = as.integer(protein_ids))
  pho2gene <- phospho$info %>%
    select(id, gene_name) %>%
    separate_rows(gene_name, sep = ";")

  # Phospho sites containing multiple genes mess up functional enrichment, as
  # these genes are usually related. As a simple hack we only take the first
  # gene on the list.
  pho2gene_first <- pho2gene %>%
    drop_na() %>%
    group_by(id) %>%
    summarise(gene_name = first(gene_name))

  list(
    pho = phospho,
    pep = peptides,
    pro = proteins,
    de = phospho_de,
    go_cc = fterms$go_cc,
    go_bp = fterms$go_bp,
    go_mf = fterms$go_mf,
    reactome = fterms$re,
    kegg = fterms$kg,
    pho2pep = pho2pep,
    pho2pro = pho2pro,
    pho2gene = pho2gene,
    pho2gene_first = pho2gene_first
  )
}
