shiny_data_de <- function(phospho, peptides, proteins, phospho_de, bm_genes, bm_go, reactome, kegg) {
  list(
    pho = phospho,
    pep = peptides,
    pro = proteins,
    de = phospho_de,
    go = bm_go,
    reactome = reactome,
    kegg = kegg
  )
}

