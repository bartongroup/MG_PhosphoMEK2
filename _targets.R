library(targets)
library(tarchetypes)

#future::plan(future.callr::callr)

packages <- c("biomaRt", "viridis", "cowplot", "ggridges", "ggbeeswarm", "GGally", "ggrepel", "uwot", "limma", "tidyverse")
tar_option_set(packages = packages, format = "qs")
options(tidyverse.quiet = TRUE, dplyr.summarise.inform = FALSE)
if(!dir.exists("tab")) dir.create("tab")

select = dplyr::select

# for interactive session only
if(interactive()) sapply(packages, library, character.only=TRUE)

files_R <- list.files(c("R", "targets"), pattern="*.R$", full.names=TRUE)
sr_ <- sapply(files_R, source)


sesinfo <- list(
  tar_target(session_info, sessionInfo())
)


c(
  sesinfo,
  targets_main()
)

