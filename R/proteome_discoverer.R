# parse site string, e.g. S477(100); S478(100)
parse_sites <- function(sites) {
  sites %>% 
    str_split(";") %>% 
    unlist() %>% 
    str_remove("\\s") %>% 
    map_dfr(function(x) {
      if (!str_detect(x, "\\(")) {
        #x <- x %>% str_remove("\\[") %>% str_remove("\\]")
        return(c(residue = x, position = NA_character_, prob = NA_character_))
      }
      # split into residue, position and probability
      v <- str_extract_all(x, "^[A-Z]|[\\d\\.]+", simplify = TRUE) %>%
        as.vector()
      names(v) <- c("residue", "position", "prob")
      v
    }) 
}


# Parse modifications string
# Returns tibble with residue, pos and prob
parse_mod <- function(s) {
  s %>% 
    str_extract("(?<=(Phospho\\s\\[))[\\w\\.\\;\\(\\)\\s\\/]+(?=(\\]))") %>%   # extract from Phospho [...]
    parse_sites() %>% 
    summarise(positions = list(position), residues = list(residue))
}

parse_prot_mod <- function(s) {
  if (is.na(s)) return(tibble(prot_positions = list(NA_real_)))
  s %>% 
    str_split("(?<=]);") %>%   # split into proteins
    unlist() %>% 
    str_remove("^\\s") %>% 
    map_dfr(function(x) {
      prot <- str_extract(x, "(^\\w+)(?=\\s)")
      str_extract(x, "(?<=(Phospho\\s\\[))[\\w\\.\\;\\(\\)\\s\\/]+(?=(\\]))") %>%    # extract from Phospho [...]  
        parse_sites() %>% 
        add_column(protein = prot) %>% 
        mutate(prot_position = str_glue("{protein}:{residue}{position}"))
    }) %>% 
    filter(!is.na(position)) %>%
    summarise(prot_positions = list(prot_position), multi = n())
}

merge_mod_info <- function(mds) {
  mds %>% 
    summarise(
      residues = str_replace_na(residue) %>% str_c(collapse = ","),
      positions = str_replace_na(position) %>% str_c(collapse = ",")
    )
}

make_modified_sequence <- function(sq, pos) {
  s <- str_to_lower(sq)
  for (i in na.omit(pos)) str_sub(s, i, i) <- str_sub(sq, i, i)
  s
}

parse_modifications <- function(d) {
  d %>% 
    rowwise() %>% 
    mutate(parse_mod(modifications)) %>% 
    mutate(parse_prot_mod(modifications_prot)) %>% 
    mutate(seqmod = make_modified_sequence(sequence, positions)) %>% 
    mutate(
      residues = str_replace_na(residues) %>% str_c(collapse = ","),
      positions = str_replace_na(positions) %>% str_c(collapse = ","),
      prot_positions = str_replace_na(prot_positions) %>% str_c(collapse = ";")
    ) %>% 
    ungroup()
}

parse_descriptions <- function(d) {
  desc <- d %>% 
    select(uniprot, description) %>% 
    distinct() %>% 
    mutate(id = row_number()) %>% 
    mutate(description = str_split(description, "\\r\\n")) %>%
    unnest(description) %>% 
    separate(description, c("prot", "rest"), " OS=") %>% 
    mutate(gene = str_extract(rest, "GN=(\\w+)") %>% str_remove("GN=")) %>% 
    select(-rest) %>% 
    group_by(id) %>% 
    summarise(description = str_c(prot, collapse = ";"), gene = str_c(gene, collapse = ";"), uniprot = first(uniprot)) %>% 
    select(-id)
  d %>% 
    select(-description) %>% 
    left_join(desc, by = "uniprot")
}

get_peptides <- function(mraw) {
  mraw %>% 
    select(
      peptide_id,
      sequence = Sequence,
      modifications = Modifications,
      modifications_prot = `Modifications in Master Proteins`,
      n_proteins = `# Proteins`,
      n_groups = `# Protein Groups`,
      pep = "Qvality PEP",
      uniprot = "Master Protein Accessions",
      description = "Master Protein Descriptions"
    )
}

normalise_peptides <- function(pep, prot, pep_info) {
  pp <- pep_info %>% 
    select(sequence, uniprot) %>% 
    distinct()
  abu <- pep$abu %>% 
    left_join(pp, by = "sequence") %>% 
    left_join(prot$abu %>% rename(prot_abu_norm = abundance_norm, prot_abu = abundance), by = c("uniprot" = "accession", "condition", "replicate")) %>% 
    drop_na() %>% 
    mutate(abundance_protnorm = mean(prot_abu_norm) * abundance_norm / prot_abu_norm)
  tab <- abu %>% 
    pivot_wider(id_cols = c("sequence", "modifications"), names_from = c("condition", "replicate"), values_from = "abundance_protnorm")
  
  pep$abu = abu
  pep$tab_protnorm = tab
  return(pep)
}

capital_mod <- function(s, pos) {
  s <- tolower(s)
  substr(s, pos, pos) <- toupper(substr(s, pos, pos))
  s
}

pad_sequences <- function(sl) {
  mx <- max(sl$position)
  ss <- rep("-", mx) %>% str_c(collapse = "")
  sl %>% 
    mutate(n_pad = max(position) - position) %>% 
    mutate(pad = str_sub(ss, 1, n_pad)) %>% 
    unite(sequence, c(pad, sequence), sep = "")
}


make_seq_list <- function(seqmod, pep_info) {
  seqmod %>% 
    select(c(sequence, modifications)) %>% 
    left_join(pep_info, by = c("sequence", "modifications")) %>% 
    filter(position > 0) %>% 
    mutate(sequence = capital_mod(sequence, position)) %>% 
    unite(label, c(uniprot, residue, position), remove = FALSE) %>% 
    mutate(label = str_remove_all(label, "\\s")) %>% 
    select(sequence, label, residue, position) %>% 
    pad_sequences()
}

write_seq_list <- function(sl, file) {
  seqinr::write.fasta(sl$sequence %>% str_split(""), sl$label, file)
}

make_pep_parse_example <- function(pepinf, seed = 123, n = 10) {
  set.seed(seed)
  seqs <- pepinf$sequence %>%
    unique() %>% 
    sample(n)
  pepinf %>% 
    filter(sequence %in% seqs) %>% 
    select(sequence, modifications, residues, positions)
}

process_pd_data <- function(raw, meta) {
  abu_rat <- raw %>% 
    select(peptide_id, contains("Abundance Ratios (log2)")) %>% 
    set_names("peptide_id", paste0("ratio_", unique(meta$replicate))) %>% 
    pivot_longer(
      cols = -peptide_id,
      names_to = c("group", "replicate"),
      values_to = "ratio",
      names_sep = "_"      
    )
  
  abu_tab <- raw %>% 
    select(peptide_id, contains("Abundance:")) %>% 
    set_names("peptide_id", meta$sample)
  
  abu_tab_norm <- raw %>% 
    select(peptide_id, contains("Abundances (Normalized)")) %>% 
    set_names("peptide_id", meta$sample)
  
  abu <- abu_tab %>% 
    pivot_longer(
      cols = meta$sample,
      names_to = c("group", "replicate"),
      values_to = "abundance",
      names_sep = "-"
    )
  
  abu_norm <- abu_tab_norm %>% 
    pivot_longer(
      cols = meta$sample,
      names_to = c("group", "replicate"),
      values_to = "abundance_norm",
      names_sep = "-"
    )
  
  
  abu <- abu %>% full_join(abu_norm, by = c("peptide_id", "group", "replicate")) %>% 
    unite("sample", group, replicate, sep = "-", remove = FALSE) %>% 
    mutate(sample = factor(sample, levels = meta$sample), group = factor(group, levels = levels(meta$group))) %>% 
    mutate(replicate = as.integer(replicate))
  
  list(abu = abu, rat = abu_rat, tab = abu_tab, tab_norm = abu_tab_norm, metadata = meta)
}

read_and_process_pd_data <- function(file, meta) {
  meta$group <- meta$condition
  raw <- readxl::read_excel(file) %>% 
    mutate(peptide_id = paste0("PEP", row_number()), .before = 1)
  peptides <- process_pd_data(raw, meta)
  peptides$info <- raw %>% 
    get_peptides() %>% 
    parse_modifications() %>% 
    parse_descriptions()
  peptides  
}


map_mq_pd <- function(mq_info, pd_info) {
  mq_pos <- mq_info %>%
    select(id, protein, amino_acid, position) %>%
    mutate(protein = str_remove(protein, "\\-\\d+")) %>% 
    mutate(prot_position = str_glue("{protein}:{amino_acid}{position}")) %>% 
    select(mq_id = id, prot_position)
  pd_pos <- pd_info %>% 
    select(pd_id = peptide_id, prot_position = prot_positions) %>% 
    separate_rows(prot_position, sep = ";")
  full_join(mq_pos, pd_pos, by = "prot_position")
}
