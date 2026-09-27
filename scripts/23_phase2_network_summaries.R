# =============================================================================
# 23_phase2_network_summaries.R · Phase 2: the numbers we report about the networks
# =============================================================================
# What we do here
#   From the fitted networks and the bootstrap we compute: the filtered density and the number of
#   significant edges with |beta| >= 0.05 (a descriptive threshold, also used to draw Figure 2); the
#   global strength of the full network; the difference in temporal global strength between arms,
#   pairing bootstrap iterations by their index; and how often each node ranks first or in the top
#   three by temporal out-strength.
#
# Input   tables/mlvar_cbwt.rds, tables/mlvar_cbt.rds
#         bootstrap files from 22_phase2_bootstrap.R
# Output  tables/network_indices.csv
#         tables/network_strength_difference.csv
#         tables/mlvar_rank_stability.csv
# Reported in  Results, Phase 2; Tables S3 and S6
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(mlVAR) })

THRESH <- 0.05
m_cbwt <- readRDS("tables/mlvar_cbwt.rds")
m_cbt <- readRDS("tables/mlvar_cbt.rds")

net_summary <- function(model, kind) {
  M_full <- getNet(model, kind, nonsig = "show")   # full network
  M_sig  <- getNet(model, kind, nonsig = "hide")   # significant edges only
  M_sig[abs(M_sig) < THRESH] <- 0
  tibble(network = kind,
         n_edges_filtered     = sum(M_sig != 0),
         filtered_density     = mean(M_sig != 0),
         global_strength_full = sum(abs(M_full)))
}

gs_boot <- read_csv("tables/mlvar_global_strength_B1000.csv", show_col_types = FALSE) %>%
  rename(gs_boot_mean = mean, gs_lo = lo, gs_hi = hi)

idx <- bind_rows(
  net_summary(m_cbwt, "temporal")        %>% mutate(arm = "CBWT"),
  net_summary(m_cbt, "temporal")        %>% mutate(arm = "CBT"),
  net_summary(m_cbwt, "contemporaneous") %>% mutate(arm = "CBWT"),
  net_summary(m_cbt, "contemporaneous") %>% mutate(arm = "CBT")
) %>% left_join(gs_boot, by = c("arm", "network")) %>%
  select(arm, network, n_edges_filtered, filtered_density, global_strength_full,
         gs_boot_mean, gs_lo, gs_hi, n_boot)
write_csv(idx, "tables/network_indices.csv")
cat("=== Network indices ===\n"); print(idx %>% mutate(across(where(is.numeric), ~ round(., 3))))

raw <- read_csv("tables/mlvar_centrality_bootstrap_B1000_raw.csv", show_col_types = FALSE)
gs <- raw %>% filter(network == "temporal") %>% distinct(arm, boot, global_strength) %>%
  pivot_wider(names_from = arm, values_from = global_strength) %>% filter(!is.na(CBWT), !is.na(CBT))
dif <- tibble(
  point_difference     = idx$global_strength_full[idx$arm == "CBWT" & idx$network == "temporal"] -
                         idx$global_strength_full[idx$arm == "CBT" & idx$network == "temporal"],
  boot_mean_difference = mean(gs$CBWT) - mean(gs$CBT),
  ci_lo = unname(quantile(gs$CBWT - gs$CBT, .025)),
  ci_hi = unname(quantile(gs$CBWT - gs$CBT, .975)),
  pct_strength_higher = 100 * mean(gs$CBWT > gs$CBT),
  n_pairs = nrow(gs))
write_csv(dif, "tables/network_strength_difference.csv")
cat("\n=== Temporal global strength, strength-based minus barrier-based ===\n")
print(dif %>% mutate(across(where(is.numeric), ~ round(., 3))))

rk <- raw %>% filter(network == "temporal") %>%
  group_by(arm, boot) %>%
  mutate(rank = rank(-out_strength, ties.method = "first")) %>%
  group_by(arm, node) %>%
  summarise(pct_first = 100 * mean(rank == 1), pct_top3 = 100 * mean(rank <= 3), .groups = "drop") %>%
  arrange(arm, desc(pct_first), desc(pct_top3))
write_csv(rk, "tables/mlvar_rank_stability.csv")
cat("\n=== Rank stability of temporal out-strength ===\n"); print(rk, n = Inf)
