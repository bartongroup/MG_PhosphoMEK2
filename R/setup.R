# Each column specification needs to contain a column "protein"
# This is used to remove rev and con

REPORTERS <- 1:10

EVIDENCE_FILE <- "mq_data/evidence.txt"

EVIDENCE_DATA_COLUMNS <- tibble::tribble(
  ~name, ~raw_name, ~type,
  "sequence", "Sequence", "c",
  "modified_sequence", "Modified sequence", "c",
  "modifications", "Modifications", "c",
  "phospho_probabilities", "Phospho (STY) Probabilities", "c",
  "proteins", "Proteins", "c",
  "protein", "Leading razor protein", "c",
  "gene_name", "Gene names", "c",
  "protein_names", "Protein names", "c",
  "experiment", "Experiment", "c",
  "charge", "Charge", "n",
  "reverse", "Reverse", "c",
  "contaminant", "Potential contaminant", "c"
)

# protein groups file

PROTEINS_FILE <- "mq_data/proteinGroups.txt"

PROTEINS_DATA_COLUMNS <- tibble::tribble(
  ~name, ~raw_name, ~type,
  "id", "id", "i",
  "peptide_ids", "Peptide IDs", "c",
  "phospho_ids", "Phospho (STY) site IDs", "c",
  "proteins", "Protein IDs", "c",
  "protein", "Majority protein IDs", "c",
  "gene_name", "Gene names", "c",
  "protein_names", "Protein names", "c",
  "sequence_length", "Sequence length", "n"
)

PROTEINS_ID_COLUMNS <- c("id")

PROTEINS_MEASURE_COLUMNS <- tibble::tibble(
  reporter = REPORTERS,
  column_name = glue::glue("Reporter intensity corrected {reporter} Phospho") |> as.character()
)

KEEP_PROTEINS_COLUMNS <- c("peptide_ids", "phospho_ids", "protein", "gene_name")



# peptides file

PEPTIDES_FILE <- "mq_data/peptides.txt"

PEPTIDES_DATA_COLUMNS <- tibble::tribble(
  ~name, ~raw_name, ~type,
  "id", "id", "i",
  "phospho_ids", "Phospho (STY) site IDs", "c",
  "protein_ids", "Protein group IDs", "c", 
  "sequence", "Sequence", "c",
  "proteins", "Proteins", "c",
  "protein", "Leading razor protein", "c",
  "gene_name", "Gene names", "c",
  "protein_names", "Protein names", "c",
  "start_position", "Start position", "n",
  "end_position", "End position", "n"
)

PEPTIDES_ID_COLUMNS <- c("id")

PEPTIDES_MEASURE_COLUMNS <- tibble::tibble(
  reporter = REPORTERS,
  column_name = glue::glue("Reporter intensity corrected {reporter}") |> as.character()
)

KEEP_PEPTIDES_COLUMNS <- c("phospho_ids", "protein", "gene_name")


# phospho file

PHOSPHO_FILE <- "mq_data/Phospho (STY)Sites.txt"

PHOSPHO_DATA_COLUMNS <- tibble::tribble(
  ~name, ~raw_name, ~type,
  "id", "id", "i",
  "peptide_ids", "Peptide IDs", "c",
  "protein_ids", "Protein group IDs", "c", 
  "proteins", "Proteins", "c",
  "protein", "Protein", "c",
  "gene_name", "Gene names", "c",
  "protein_names", "Protein names", "c",
  "localization_prob", "Localization prob", "n",
  "score_for_locatization", "Score for localization", "n",
  "num_phospho", "Number of Phospho (STY)", "n",
  "amino_acid", "Amino acid", "c",
  "sequence_window", "Sequence window", "c",
  "modification_window", "Modification window", "c",
  "peptide_window_coverage", "Peptide window coverage", "c",
  "phospho_probabilities", "Phospho (STY) Probabilities", "c",
  "position_in_peptide", "Position in peptide", "n",
  "position", "Position", "n",
  "charge", "Charge", "n"
)

PHOSPHO_ID_COLUMNS <- c("id")

KEEP_PHOSPHO_COLUMNS <- c("protein", "gene_name", "position", "localization_prob")

PHOSPHO_MEASURE_COLUMNS <- tibble::tibble(
  reporter = REPORTERS,
  column_name = glue::glue("Reporter intensity corrected {reporter}___1") |> as.character()
)


SAMPLE_REPORTER <- tibble::tibble(
  reporter = REPORTERS,
  sample = c("DMSO-1", "DMSO-2", "DMSO-3", "DMSO-4", "DMSO-5", "MEKi-1", "MEKi-2", "MEKi-3", "MEKi-4", "MEKi-5")
)

CONDITIONS <- tibble::tibble(CONDITION = c("DMSO", "MEKi"))

NORMALISATIONS <- tibble::tibble(
  NAME = c("median", "constand", "protein", "mean_protein"),
  WHAT = c("value_med", "value_constand", "value_prot", "value_prot_mean"),
  LOG = c(TRUE, FALSE, TRUE, TRUE)
)

NORM_COND <- tidyr::expand_grid(NORMALISATIONS, CONDITIONS) |>
  dplyr::mutate(NAME = paste0(NAME, "_", CONDITION))

FDR_LIMIT <- 0.05
LOGFC_LIMIT <- 0
