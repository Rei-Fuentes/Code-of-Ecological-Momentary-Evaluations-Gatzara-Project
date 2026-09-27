# =============================================================================
# 43_phase4_cluster_stability.R · Phase 4: how stable are the clusters?
# =============================================================================
# What we do here
#   We resample participants with replacement within each arm 30 times, refit k-means with the chosen
#   number of clusters and record the mean silhouette. Duplicated participants sit at zero distance
#   from each other, which inflates the silhouette, so these values are not evidence of better
#   separation than in the original sample.
#
# Input   tables/idionomic_betas.csv
#         tables/idionomic_clusters.csv
# Output  tables/cluster_bootstrap_silhouette.csv
# Reported in  Supplementary S8
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(cluster) })

set.seed(20260519)

idio  <- read_csv("tables/idionomic_betas.csv", show_col_types = FALSE)
clust <- read_csv("tables/idionomic_clusters.csv", show_col_types = FALSE)

mat <- idio %>%
  mutate(pair = paste(outcome, process, sep = "__")) %>%
  select(code, arm, pair, beta) %>%
  pivot_wider(names_from = pair, values_from = beta) %>%
  filter(if_all(-c(code, arm), ~ !is.na(.)))

boot_silhouette <- function(arm_label, n_boot = 30) {
  sub <- mat %>% filter(arm == arm_label)
  k_best <- clust %>% filter(arm == arm_label) %>% pull(silhouette_best_k) %>% unique()
  replicate(n_boot, {
    idx <- sample(seq_len(nrow(sub)), replace = TRUE)
    Xb <- scale(as.matrix(sub[idx, -(1:2)]))
    Xb <- Xb[apply(Xb, 1, function(r) all(is.finite(r))), , drop = FALSE]
    if (nrow(Xb) < k_best + 1) return(NA)
    clb <- kmeans(Xb, centers = k_best, nstart = 10)$cluster
    mean(silhouette(clb, dist(Xb))[, 3])
  })
}

res <- bind_rows(lapply(c("CBWT", "CBT"), function(a) {
  s <- boot_silhouette(a)
  tibble(arm = a, n_participants = sum(mat$arm == a), n_boot = length(s),
         mean_silhouette = mean(s, na.rm = TRUE), sd = sd(s, na.rm = TRUE),
         min = min(s, na.rm = TRUE), max = max(s, na.rm = TRUE))
}))
write_csv(res, "tables/cluster_bootstrap_silhouette.csv")
print(res %>% mutate(across(where(is.numeric), ~ round(., 3))))
