# =============================================================================
# 04_psychometrics.R · How reliable are the single daily items?
# =============================================================================
# What we do here
#   Each construct is measured with one item, so classic internal consistency does not apply. We
#   report the indices recommended for EMA (Cranford et al., 2006; Bolger & Laurenceau, 2013) on the
#   intervention-phase prompts: the share of variance between people (ICC1), the reliability of a
#   participant's mean (Spearman-Brown with the median number of prompts), the lag-1 autocorrelation
#   within people, and the two-item coefficient of the global well-being composite.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/psychometrics_icc.csv
#         tables/psychometrics_autocorrelations.csv
#         tables/psychometrics_composite.csv
# Reported in  Method, EMA item development and reliability; Supplementary S2
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(lme4)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE) %>%
  filter(phase == "intervention")

items <- c("hedonic_wellbeing","cognitive_wellbeing","emotional_awareness",
           "compassion","self_compassion","positive_emotion_regulation",
           "gratitude","negative_emotion_raw","cognitive_reappraisal")

# --- ICC1 and R_kr per item ---
icc_results <- tibble()
for (it in items) {
  d <- df %>% select(code, value = all_of(it)) %>% filter(!is.na(value))
  fit <- lmer(value ~ 1 + (1|code), data = d, REML = TRUE)
  vc <- as.data.frame(VarCorr(fit))
  tau2 <- vc$vcov[vc$grp == "code"]
  sigma2 <- vc$vcov[vc$grp == "Residual"]
  icc1 <- tau2 / (tau2 + sigma2)
  # median observations per person (for the Spearman-Brown R_kr)
  n_per_person <- d %>% group_by(code) %>% summarise(n = n()) %>% pull(n)
  k <- median(n_per_person)
  r_kr <- (k * icc1) / (1 + (k - 1) * icc1)
  icc_results <- bind_rows(icc_results, tibble(
    item = it, n_obs = nrow(d), n_subj = n_distinct(d$code),
    median_obs_per_person = k,
    ICC1 = icc1, R_kr = r_kr
  ))
}
print(icc_results %>% mutate(across(c(ICC1, R_kr), ~round(.x, 3))))

# --- Two-item coefficient for global_wellbeing (hedonic + cognitive) ---
gw <- df %>% select(hedonic_wellbeing, cognitive_wellbeing) %>%
  drop_na() %>% as.matrix()
# Spearman-Brown 2-item: 2*r / (1+r)
r12 <- cor(gw)[1,2]
alpha_2item <- (2 * r12) / (1 + r12)
cat(sprintf("\n=== global_wellbeing composite (2 items) ===\n"))
cat(sprintf("Hedonic-cognitive correlation across prompts: r = %.3f\n", r12))
cat(sprintf("Spearman-Brown 2-item alpha:                %.3f\n", alpha_2item))
write_csv(tibble(inter_item_r = r12, spearman_brown = alpha_2item), "tables/psychometrics_composite.csv")

# --- Lag-1 autocorrelation per person and item ---
cat("\n=== Lag-1 within-person autocorrelation per item (median and IQR) ===\n")
ac_results <- tibble()
for (it in items) {
  per_person_ac <- df %>%
    arrange(code, date) %>%
    group_by(code) %>%
    summarise(ac = ifelse(sum(!is.na(.data[[it]])) >= 5,
                          cor(.data[[it]], dplyr::lag(.data[[it]]),
                              use = "pairwise.complete.obs"),
                          NA_real_),
              .groups = "drop") %>%
    filter(!is.na(ac))
  ac_results <- bind_rows(ac_results, tibble(
    item = it,
    n_subj = nrow(per_person_ac),
    ac_median = median(per_person_ac$ac, na.rm = TRUE),
    ac_q25 = quantile(per_person_ac$ac, .25, na.rm = TRUE),
    ac_q75 = quantile(per_person_ac$ac, .75, na.rm = TRUE)
  ))
}
print(ac_results %>% mutate(across(starts_with("ac_"), ~round(.x, 3))))

# Save
write_csv(icc_results, "tables/psychometrics_icc.csv")
write_csv(ac_results,  "tables/psychometrics_autocorrelations.csv")
cat("\nWritten: tables/psychometrics_icc.csv\n")
cat("Written: tables/psychometrics_autocorrelations.csv\n")
