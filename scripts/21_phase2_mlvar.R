# =============================================================================
# 21_phase2_mlvar.R · Phase 2: daily networks in each arm
# =============================================================================
# What we do here
#   We estimate a multilevel vector autoregression separately in each arm (Sahdra et al., 2025), with
#   the nine daily items and the participants who have at least 20 intervention-phase prompts. It
#   gives a temporal network (how today's state predicts tomorrow's), a contemporaneous network (how
#   states go together on the same day) and a between-person network.
#   Settings: lags = 1, estimator = "lmer", orthogonal temporal and contemporaneous effects. We pass
#   no day variable, so each lag links a participant's successive prompts, including those separated
#   by more than one day (12.8% of lags).
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/mlvar_cbwt.rds, tables/mlvar_cbt.rds
#         tables/mlvar_centrality.csv, tables/mlvar_edges.csv, tables/mlvar_summary.csv
# Reported in  Results, Phase 2; Figure 2; Supplementary S6
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(mlVAR); library(qgraph)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
d  <- df %>% filter(phase == "intervention")

vars <- c("hedonic_wellbeing","cognitive_wellbeing","emotional_awareness",
          "compassion","self_compassion","positive_emotion_regulation",
          "gratitude","negative_emotion_regulation","cognitive_reappraisal")

# beep and day are computed for reference only; they are not passed to mlVAR.
d <- d %>%
  arrange(code, date) %>%
  group_by(code) %>%
  mutate(beep = row_number(), day = as.integer(date - min(date))) %>%
  ungroup()

cbwt <- d %>% filter(arm == "CBWT")
cbt <- d %>% filter(arm == "CBT")

# Participants with >= 20 intervention-phase observations.
keep <- function(x) x %>% group_by(code) %>% filter(n() >= 20) %>% ungroup()
cbwt <- keep(cbwt); cbt <- keep(cbt)

cat(sprintf("CBWT: %d obs / %d participants\n", nrow(cbwt), n_distinct(cbwt$code)))
cat(sprintf("CBT: %d obs / %d participants\n", nrow(cbt), n_distinct(cbt$code)))

fit_mlvar <- function(dat) {
  # Sequential order, no dayvar (see header).
  mlVAR(
    data = as.data.frame(dat),
    vars = vars,
    idvar = "code",
    lags = 1,
    estimator = "lmer",
    temporal = "orthogonal",
    contemporaneous = "orthogonal",
    verbose = FALSE
  )
}

cat("\nFitting mlVAR, strength-based arm (CBWT)...\n")
m_cbwt <- fit_mlvar(cbwt)
cat("Fitting mlVAR, barrier-based arm (CBT)...\n")
m_cbt <- fit_mlvar(cbt)

saveRDS(m_cbwt, "tables/mlvar_cbwt.rds")
saveRDS(m_cbt, "tables/mlvar_cbt.rds")

# Networks and centrality (in-strength + out-strength = strength).
get_mats <- function(m) {
  list(
    temporal      = getNet(m, "temporal", nonsig = "show"),
    contemp       = getNet(m, "contemporaneous", rule = "and", nonsig = "show"),
    between       = getNet(m, "between", rule = "and", nonsig = "show")
  )
}
mats_cbwt <- get_mats(m_cbwt)
mats_cbt <- get_mats(m_cbt)

cent <- function(M, lab) {
  cm <- centrality(M)
  tibble(node = vars, type = lab,
         in_strength  = cm$InDegree,
         out_strength = cm$OutDegree,
         strength     = cm$OutDegree + cm$InDegree)
}

cent_all <- bind_rows(
  cent(mats_cbwt$temporal, "CBWT_temporal"),
  cent(mats_cbt$temporal, "CBT_temporal"),
  cent(mats_cbwt$contemp,  "CBWT_contemp"),
  cent(mats_cbt$contemp,  "CBT_contemp")
)
write_csv(cent_all, "tables/mlvar_centrality.csv")

# Full matrices in long format.
to_long <- function(M, lab) {
  as.data.frame(M) %>%
    mutate(from = rownames(M)) %>%
    pivot_longer(-from, names_to = "to", values_to = "weight") %>%
    mutate(network = lab)
}
edges <- bind_rows(
  to_long(mats_cbwt$temporal, "CBWT_temporal"),
  to_long(mats_cbt$temporal, "CBT_temporal"),
  to_long(mats_cbwt$contemp,  "CBWT_contemp"),
  to_long(mats_cbt$contemp,  "CBT_contemp")
)
write_csv(edges, "tables/mlvar_edges.csv")

# Unfiltered density and global strength (sum of absolute edge weights, full network).
density <- function(M) mean(M != 0)
strength_global <- function(M) sum(abs(M))

summ <- tibble(
  network = c("CBWT_temporal","CBT_temporal","CBWT_contemp","CBT_contemp"),
  density = c(density(mats_cbwt$temporal), density(mats_cbt$temporal),
              density(mats_cbwt$contemp),  density(mats_cbt$contemp)),
  global_strength = c(strength_global(mats_cbwt$temporal), strength_global(mats_cbt$temporal),
                      strength_global(mats_cbwt$contemp),  strength_global(mats_cbt$contemp))
)
write_csv(summ, "tables/mlvar_summary.csv")

cat("\n=== mlVAR resumen ===\n")
print(summ)
cat("\nTop three nodes by strength per network:\n")
print(cent_all %>% group_by(type) %>% slice_max(strength, n = 3))
cat("\nWritten to tables/.\n")
