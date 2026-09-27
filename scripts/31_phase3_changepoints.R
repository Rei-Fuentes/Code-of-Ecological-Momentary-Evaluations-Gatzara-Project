# =============================================================================
# 31_phase3_changepoints.R · Phase 3: when do people change?
# =============================================================================
# What we do here
#   Following Snippe et al. (2024), we look for sudden shifts in each participant's daily series with
#   the E-divisive method (ecp): alpha = 2, min.size = 14, 1,000 permutations, significance .005
#   (Bonferroni for ten variables), in series with at least 28 prompts (2 x min.size).
#   A change-point counts as a gain when the mean after it is higher than before it and the whole
#   series is improving, judged by a second pass with a single change-point (k = 1, min.size =
#   max(7, n/3), 199 permutations, significance .05). Cohen's d compares the 14 prompts after and
#   the 14 before each change-point. Finally we test, with a one-sided binomial test per arm, whether
#   a process tends to improve before an outcome (more than 7 prompts apart), with the FDR within arm.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/changepoints_all.csv
#         tables/changepoints_summary.csv
#         tables/changepoint_sequence_tests.csv
# Reported in  Results, Phase 3; Table S4; Figure S3; Supplementary S7
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ecp)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
d  <- df %>% filter(phase == "intervention") %>% arrange(code, date)

vars <- c("hedonic_wellbeing","cognitive_wellbeing","global_wellbeing","emotional_awareness",
          "compassion","self_compassion","positive_emotion_regulation",
          "gratitude","negative_emotion_regulation","cognitive_reappraisal")


run_edivisive <- function(y, min_size = 14, R = 1000, sig = .005) {
  if (length(y) < min_size * 2) return(NULL)
  fit <- e.divisive(matrix(y, ncol = 1), sig.lvl = sig, R = R, k = NULL,
                    min.size = min_size, alpha = 2)
  cps <- fit$estimates
  # estimates include the series ends; keep interior change-points only.
  cps <- cps[cps > 1 & cps < length(y)]
  list(cps = cps, n_cluster = fit$k.hat)
}

# Overall improvement: single change-point pass (k = 1).
trajectory_improves <- function(y) {
  if (length(y) < 14) return(FALSE)
  fit <- tryCatch(
    e.divisive(matrix(y, ncol = 1), sig.lvl = .05, R = 199, k = 1,
               min.size = max(7, floor(length(y)/3)), alpha = 2),
    error = function(e) NULL
  )
  if (is.null(fit) || length(fit$estimates) < 3) return(FALSE)
  cp <- fit$estimates[2]
  if (is.na(cp) || cp < 2 || cp >= length(y)) return(FALSE)
  mean(y[cp:length(y)], na.rm = TRUE) > mean(y[1:(cp-1)], na.rm = TRUE)
}

# Change-points per participant and variable, each classified as gain or not.
results <- list(); i <- 1
for (id in unique(d$code)) {
  sub <- d %>% filter(code == id)
  arm <- sub$arm[1]
  for (v in vars) {
    y <- sub[[v]]
    if (sum(!is.na(y)) < 28) next   # 2 x min.size
    out <- tryCatch(run_edivisive(y), error = function(e) NULL)
    if (is.null(out) || length(out$cps) == 0) next
    overall_up <- trajectory_improves(y)
    for (cp in out$cps) {
      w_before <- y[max(1, cp - 14):(cp - 1)]
      w_after  <- y[cp:min(length(y), cp + 13)]
      delta    <- mean(w_after, na.rm = TRUE) - mean(w_before, na.rm = TRUE)
      sd_pool  <- sqrt((sd(w_before, na.rm = TRUE)^2 + sd(w_after, na.rm = TRUE)^2) / 2)
      d_eff    <- if (is.finite(sd_pool) && sd_pool > 0) delta / sd_pool else NA_real_
      is_gain  <- !is.na(delta) && delta > 0 && overall_up
      results[[i]] <- tibble(code = id, arm = arm, variable = v,
                             cp_index = cp, delta = delta, cohen_d = d_eff,
                             overall_improving = overall_up, is_gain = is_gain)
      i <- i + 1
    }
  }
}
cps <- bind_rows(results)
write_csv(cps, "tables/changepoints_all.csv")

# Summary by variable and arm.
gains <- cps %>% filter(is_gain)
summary_gains <- gains %>%
  group_by(arm, variable) %>%
  summarise(n_gains = n(),
            n_subj_with_gain = n_distinct(code),
            median_d = median(cohen_d, na.rm = TRUE),
            .groups = "drop")
write_csv(summary_gains, "tables/changepoints_summary.csv")

cat("\n=== Gains by arm ===\n")
print(summary_gains)

# Binomial test: does the first gain in a process precede the first gain in an outcome?
outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")

first_gain <- gains %>%
  group_by(code, arm, variable) %>%
  summarise(first_cp = min(cp_index), .groups = "drop")

binom_rows <- list(); j <- 1
for (br in c("CBWT", "CBT")) {
  for (p in processes) {
    for (o in outcomes) {
      sub <- first_gain %>% filter(arm == br, variable %in% c(p, o))
      pp  <- sub %>% filter(variable == p) %>% select(code, cp_p = first_cp)
      oo  <- sub %>% filter(variable == o) %>% select(code, cp_o = first_cp)
      pair <- inner_join(pp, oo, by = "code")
      # Precedence requires a separation of more than 7 observations
      pair <- pair %>% mutate(diff = cp_p - cp_o,
                              po = diff < -7,   # process before outcome
                              op = diff >  7)   # outcome before process
      n_po <- sum(pair$po); n_op <- sum(pair$op)
      n    <- n_po + n_op
      if (n == 0) next
      bt <- binom.test(n_po, n, p = .5, alternative = "greater")
      binom_rows[[j]] <- tibble(arm = br, process = p, outcome = o,
                                n_pairs = nrow(pair), n_po = n_po, n_op = n_op,
                                p_value = bt$p.value)
      j <- j + 1
    }
  }
}
binom <- bind_rows(binom_rows)
if (nrow(binom) > 0) {
  binom <- binom %>% group_by(arm) %>%
    mutate(p_fdr = p.adjust(p_value, method = "BH")) %>% ungroup()
}
write_csv(binom, "tables/changepoint_sequence_tests.csv")

cat("\n=== Binomial precedence tests (process before outcome) ===\n")
if (nrow(binom) == 0) {
  cat("No co-occurring pairs.\n")
} else {
  print(binom %>% filter(p_fdr < .1) %>% arrange(arm, p_fdr))
}
cat("\nWritten to tables/.\n")
