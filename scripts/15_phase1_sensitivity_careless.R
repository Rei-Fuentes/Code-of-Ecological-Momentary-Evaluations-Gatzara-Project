# =============================================================================
# 15_phase1_sensitivity_careless.R · Phase 1 without careless responses
# =============================================================================
# What we do here
#   We drop the prompts in which at least eight of the nine items had the same value (longstring >= 8,
#   from 05_careless_responding.py), rebuild the triads and refit Phase 1.
#
# Input   data/processed/ema_analysis_ready.csv
#         tables/careless_responses_per_prompt.csv
# Output  tables/sensitivity_careless.csv
# Reported in  Results, Phase 1; Table S9
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(lme4); library(lmerTest); library(MuMIn)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
cr <- read_csv("tables/careless_responses_per_prompt.csv", show_col_types = FALSE)

cat(sprintf("Prompts: %d\n", nrow(df)))
cat(sprintf("Prompts flagged (longstring>=8): %d\n", sum(cr$flag_longstring)))

# Flagged prompts, identified by (code, date)
flagged <- cr %>% filter(flag_longstring == 1) %>% select(code, date)
df_clean <- df %>% anti_join(flagged, by = c("code","date"))
cat(sprintf("Prompts kept after removing flagged prompts: %d\n", nrow(df_clean)))

outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")
vars <- unique(c(outcomes, processes))

baseline <- df %>% filter(phase == "baseline") %>%
  group_by(code) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_base"),
            .groups = "drop")
cwc <- function(x) x - mean(x, na.rm = TRUE)

run_models <- function(data) {
  d <- data %>% filter(phase == "intervention") %>%
    left_join(baseline, by = "code") %>%
    arrange(code, date) %>% group_by(code) %>%
    mutate(across(all_of(vars), cwc, .names = "{.col}_c")) %>% ungroup()
  results <- list(); i <- 1
  for (Y in outcomes) for (X in processes) {
    dat <- d %>%
      arrange(code, date) %>% group_by(code) %>%
      mutate(
        Y_t   = !!sym(paste0(Y,"_c")),
        Y_tp1 = lead(!!sym(paste0(Y,"_c")), 1),
        X_t   = !!sym(paste0(X,"_c")),
        X_tm1 = lag(!!sym(paste0(X,"_c")), 1),
        Y_base = !!sym(paste0(Y,"_base")),
        d_prev = as.integer(date - lag(date)),
        d_next = as.integer(lead(date) - date)
      ) %>% ungroup() %>%
      filter(d_prev == 1, d_next == 1) %>%
      filter(!is.na(Y_t), !is.na(Y_tp1), !is.na(X_t), !is.na(X_tm1), !is.na(Y_base)) %>%
      mutate(arm = factor(arm, levels = c("CBT","CBWT")))
    if (nrow(dat) < 50) next
    fit <- tryCatch(
      lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + Y_base + (1 | code), data = dat, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$coefficients
    pull <- function(row) if (row %in% rownames(s)) s[row, c("Estimate","Std. Error","Pr(>|t|)")] else c(NA,NA,NA)
    bi <- pull("X_tm1:armCBWT")
    results[[i]] <- tibble(outcome=Y, process=X, n_obs=nrow(dat),
                           b_interact=bi[1], p_interact=bi[3])
    i <- i + 1
  }
  bind_rows(results) %>% mutate(p_fdr = p.adjust(p_interact, method="BH"))
}

cat("\n--- All prompts ---\n")
res_full <- run_models(df)
cat("\n--- Without prompts with longstring >= 8 ---\n")
res_sens <- run_models(df_clean)

cmp <- res_full %>%
  select(outcome, process, b_full=b_interact, pfdr_full=p_fdr) %>%
  left_join(res_sens %>% select(outcome, process,
                                b_sens=b_interact, pfdr_sens=p_fdr),
            by = c("outcome","process")) %>%
  mutate(db = b_sens - b_full,
         sign_flip = sign(b_full) != sign(b_sens),
         sig_full = pfdr_full < .05, sig_sens = pfdr_sens < .05)

dir.create("tables", showWarnings = FALSE)
write_csv(cmp, "tables/sensitivity_careless.csv")

cat("\n=== All prompts vs without flagged prompts ===\n")
cat(sprintf("Interactions with FDR p < .05, all=%d | without flagged=%d\n",
            sum(cmp$sig_full,na.rm=T), sum(cmp$sig_sens,na.rm=T)))
cat(sprintf("Sign changes: %d\n", sum(cmp$sign_flip,na.rm=T)))
cat(sprintf("Correlation of coefficients: %.3f\n", cor(cmp$b_full, cmp$b_sens, use="complete.obs")))
cat(sprintf("Median absolute change: %.4f | maximum: %.4f\n",
            median(abs(cmp$db),na.rm=T), max(abs(cmp$db),na.rm=T)))
cat("\nWritten: tables/sensitivity_careless.csv\n")
