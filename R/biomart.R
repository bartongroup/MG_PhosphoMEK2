biomart_fetch_genes <- function(mart) {
  getBM(attributes = c(
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
    "percentage_gene_gc_content",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id",
    "description"
  ), mart = mart) %>%
    dplyr::rename(
      chr = chromosome_name,
      start = start_position,
      end = end_position,
      gene_id = ensembl_gene_id,
      gene_name = external_gene_name,
      ncbi_id = entrezgene_id,
      gc_content = percentage_gene_gc_content
    ) %>%
    dplyr::mutate(description = str_remove(description, "\\s\\[.*\\]")) %>%
    tibble::as_tibble()
}



# For a given list of gene names, fetch all GO-terms

bm_fetch_go_genes <- function(mart, gene_names, slim = FALSE) {
  id <- ifelse(slim, "goslim_goa_accession", "go_id")
  gene2go <- getBM(
    attributes = c("external_gene_name", id),
    filters = "external_gene_name",
    values = gene_names,
    mart = mart
  ) %>%
    dplyr::rename(gene_name = external_gene_name, term_id = !!sym(id)) %>%
    dplyr::filter(term_id != "") %>%
    tibble::as_tibble()
}

# Get all GO-term descriptions

bm_fetch_go_descriptions <- function(mart) {
  # filtering on GO-terms does not work properly, so I have to fetch all terms
  getBM(
    attributes = c("go_id", "name_1006", "namespace_1003"),
    mart = mart) %>%
    dplyr::rename(
      term_id = go_id,
      term_name = name_1006,
      #term_description = definition_1006,
      term_domain = namespace_1003
    ) %>%
    dplyr::filter(term_id != "") %>%
    tibble::as_tibble()
}

# Get ontology directly from geneontology.org

go_fetch_go_descriptions <- function(obo_file = "go.obo") {
  if (!dir.exists("cache")) dir.create("cache")
  obo_path <- file.path("cache", obo_file)
  if (!file.exists(obo_path)) download.file("http://purl.obolibrary.org/obo/go.obo", obo_path)
  go <- ontologyIndex::get_ontology(obo_path, extract_tags = c("everything"))
  tibble(
    term_id = go$id,
    term_name = go$name,
    term_namespace = unlist(go$namespace)
  )
}

# Fetch GO-gene and GO descriptions

bm_fetch_go <- function(mart, gene_names, slim = FALSE) {
  gene2go <- bm_fetch_go_genes(mart, gene_names, slim)
  # using geneontology.org as Ensembl has term descriptions missing
  # goterms <- go_fetch_go_descriptions()
  goterms <- bm_fetch_go_descriptions(mart)
  terms <- gene2go$term_id %>%
    unique()
  go2gene <- map(terms, function(trm) gene2go[gene2go$term_id == trm, ]$gene_name) %>%
    set_names(terms)

  list(
    term2gene = go2gene,
    gene2term = gene2go,
    terms = goterms %>% filter(term_id %in% terms)  # full OBO table is huge
  )
}




# Reactome

reactome_fetch_genes <- function(mart, gene_names) {
  gene2re <- getBM(
    attributes = c("external_gene_name", "reactome"),
    filters = "external_gene_name",
    values = gene_names,
    mart = mart
  ) %>%
    dplyr::rename(gene_name = external_gene_name, term_id = "reactome") %>%
    dplyr::filter(term_id != "") %>%
    tibble::as_tibble()
}


reactome_fetch_pathways <- function() {
  url <- "https://reactome.org/download/current/ReactomePathways.txt"
  colms <- c("reactome_id", "name", "species")
  read_tsv(url, col_names = colms, col_types = cols())
}



# Get Reactome data in the same format as GO-data

fetch_reactome <- function(mart, gene_names) {
  r <- reactome_fetch_pathways()
  g2r <- reactome_fetch_genes(mart, gene_names)
  terms <- g2r$term_id %>%
    unique()

  reactometerms <- select(r, reactome_id, name) %>%
    rename(term_id = reactome_id, term_name = name) %>%
    filter(term_id %in% terms) %>%
    distinct()
  r2g <- map(terms, function(trm) g2r[g2r$term_id == trm, ]$gene_name) %>%
    set_names(terms)
  list(
    gene2term = g2r,
    term2gene = r2g,
    terms = reactometerms
  )
}
