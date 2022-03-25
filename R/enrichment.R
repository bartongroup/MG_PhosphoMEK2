# Enrichment function
#
# Performs hypergeometric test.
# Input: genes_sel and genes_all are vectors with gene IDs.
# gene2term: data frame with columns id and term
# gene2name: named vector to change genes ids into gene names
# term_info: data frame with columns term and other columns with name, description and so on.
# These additional columns will be inluded in the output.

functionalEnrichment <- function(genes_all, genes_sel, term_data, gene2name = NULL,
                                 min_count = 3, sig_limit = 0.05) {

  gene2term <- term_data$gene2term
  term_info <- term_data$terms

  # select only terms represented in our gene set
  gene2term <- gene2term %>%
    filter(gene_name %in% genes_all)

  # all terms present in the selection
  terms <- gene2term %>%
    filter(gene_name %in% genes_sel) %>%
    pull(term_id) %>%
    unique()

  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)

  # empty line for missing terms
  na_term <- term_info %>% slice(1) %>% mutate_all(~NA)

  res <- map_dfr(terms, function(term) {
    info <- term_info %>%
      filter(term_id == term)
    # returns NAs if no term found
    if (nrow(info) == 0) info <- na_term %>% mutate(term_id = term)

    # all genes with the term
    tgenes <- gene2term %>%
      filter(term_id == term) %>%
      pull(gene_name)
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)

    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)

    expected <- nuni * Nsel / Nuni
    fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow = 2)
    ft <- fisher.test(fish, alternative = "greater")
    p <- as.numeric(ft$p.value)

    if (!is.null(gene2name)) tgenes_sel <- gene2name[tgenes_sel] %>% unname()

    bind_cols(
      info,
      tibble(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        ids = paste(tgenes_sel, collapse = ","),
        P = p
      )
    )
  }) %>%
    mutate(P = p.adjust(P, method = "BH")) %>%
    filter(sel >= min_count & P <= sig_limit) %>%
    arrange(desc(enrich)) %>%
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))

  res
}

make_term_list <- function(gene2term) {
  gene2term %>%
    arrange(term_id) %>%
    group_split(term_id) %>%
    map(function(w) w$gene_name) %>%
    set_names(unique(gene2term$term_id) %>% sort)  # dodgy!
}

fgsea_run <- function(trm, res, min.size = 3) {
  res <- res %>%
    filter(!is.na(value) & !is.na(gene_name))
  term_list <-  trm$gene2term %>%
    filter(gene_name %in% res$gene_name) %>%
    make_term_list()
  ranks <-  set_names(res$value, res$gene_name)
  fgsea::fgsea(pathways = term_list, stats = ranks, nproc = 6, minSize = min.size, eps = 0) %>%
    as_tibble %>%
    left_join(trm$terms, by = c("pathway" = "term_id")) %>%
    arrange(NES) %>%
    select(term = pathway, term_name, pval, padj, NES, size, leading_edge = leadingEdge)
}

fgsea_cache <- function(d, terms, file, valvar = "logFC", groupvar = "contrast") {
  if (file.exists(file)) {
    fg <- read_rds(file)
  } else {
    fg <- d %>%
      mutate(value = !!sym(valvar)) %>%
      group_split(!!sym(groupvar)) %>%
      map_dfr(function(w) {
        fgsea_run(terms, w) %>%
          mutate(!!groupvar := first(w[[groupvar]]))
      })
    write_rds(fg, file)
  }
  fg
}

fgsea_all_terms <- function(d, all_terms, valvar = "logFC", groupvar = "contrast") {
  nms <- names(all_terms)
  map(nms, function(trm) {
    cat(str_glue("  Computing fgsea for {trm}\n\n"))
    cache_file <- file.path("cache", str_glue("fgsea_{trm}.rds"))
    fgsea_cache(d, all_terms[[trm]], cache_file, valvar, groupvar)
  }) %>%
    set_names(nms)
}

plot_fgsea_enrichment <- function(term, res, terms, value = "logFC") {
  lst <- terms$term2gene[[term]]
  rnks <- set_names(res[[value]], res$gene_name)
  fgsea::plotEnrichment(lst, rnks)
}

select_star_fgsea <- function(se, fg, groupvar = "contrast") {
  fg %>%
    filter(padj < 0.05) %>%
    group_split(term, !!sym(groupvar)) %>%
    map_dfr(function(w) {
      term <- as.character(w$term)
      gr <- as.character(w[[groupvar]])
      genes <- w$leading_edge[[1]]
      se %>%
        filter(gene_name %in% genes & !!sym(groupvar) == gr) %>%
        add_column(term_id = term, .before = "gene_name") %>%
        add_column(NES = w$NES, .before = "gene_name")
    })
}

select_star_go <- function(se, bm_go, terms) {
  bm_go$gene2term %>%
    filter(term_id %in% terms) %>%
    inner_join(se, by = "gene_name")
}
