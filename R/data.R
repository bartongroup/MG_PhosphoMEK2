
make_metadata <- function(samples, conditions) {
  samples %>%
    separate(sample, c("condition", "replicate"), sep = "-", remove = FALSE) %>%
    mutate(
      condition = factor(condition, levels = conditions),
      replicate = factor(replicate)
    ) %>%
    arrange(condition, replicate)
}

dat2mat <- function(d, what) {
  d %>%
    select(id, multi, sample, val = !!what) %>%
    unite(mid, c(id, multi), sep = "-") %>%
    drop_na() %>%
    pivot_wider(id_cols = mid, names_from = sample, values_from = val) %>%
    column_to_rownames("mid") %>%
    as.matrix()
}

