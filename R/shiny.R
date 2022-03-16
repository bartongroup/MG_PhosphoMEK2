shiny_data_de <- function(phospho, peptides, proteins, phospho_de, bm_genes, terms) {
  pho2pep <- phospho$info %>% 
    select(id, peptide_ids) %>% 
    separate_rows(peptide_ids, sep = ";")
  pho2pro <- phospho$info %>% 
    select(id, protein_ids) %>% 
    separate_rows(protein_ids, sep = ";")
  pho2gene <- phospho$info %>% 
    select(id, gene_name) %>% 
    separate_rows(gene_name, sep = ";")
  
  list(
    pho = phospho,
    pep = peptides,
    pro = proteins,
    de = phospho_de,
    go = terms$go,
    reactome = terms$re,
    kegg = terms$kg,
    pho2pep = pho2pep,
    pho2pro = pho2pro,
    pho2gene = pho2gene
  )
}

