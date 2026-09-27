# =============================================================================
# 13_phase1_sensitivity_min30obs.R · Phase 1 without short series
# =============================================================================
# What we do here
#   The main models keep participants with few prompts and let partial pooling shrink them. Here we
#   check that they do not drive the result: we refit Phase 1 only with participants who have at
#   least 30 intervention-phase prompts and compare the 24 interaction coefficients.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/sensitivity_min30obs.csv
# Reported in  Results, Phase 1; Table S7
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(lme4); library(lmerTest); library(MuMIn)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)

# --- Participants with >= 30 intervention-phase observations ---
n_dur <- df %>% filter(phase == "intervention") %>% count(code, name = "n_dur")
keep_ids <- n_dur %>% filter(n_dur >= 30) %>% pull(code)
drop_ids <- n_dur %>% filter(n_dur <  30) %>% pull(code)
cat(sprintf("Participants with intervention-phase data: %d\n", nrow(n_dur)))
cat(sprintf("Kept (>= 30 obs): %d | Dropped (< 30 obs): %d\n",
            length(keep_ids), length(drop_ids)))

outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")
vars <- unique(c(outcomes, processes))

# Baseline level per participant (mean in the baseline phase).
baseline <- df %>% filter(phase == "baseline") %>%
  group_by(code) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_base"),
            .groups = "drop")

cwc <- function(x) x - mean(x, na.rm = TRUE)

run_models <- function(ids) {
  d <- df %>% filter(phase == "intervention", code %in% ids) %>%
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
    results[[i]] <- tibble(
      outcome = Y, process = X, n_obs = nrow(dat), n_subj = n_distinct(dat$code),
      b_interact = bi[1], se_interact = bi[2], p_interact = bi[3]
    )
    i <- i + 1
  }
  bind_rows(results) %>%
    mutate(p_interact_fdr = p.adjust(p_interact, method = "BH")) %>%
    arrange(outcome, process)
}

cat("\n--- Full sample ---\n")
res_full <- run_models(n_dur$code)
cat("\n--- Restricted: >= 30 intervention-phase observations ---\n")
res_sens <- run_models(keep_ids)

cmp <- res_full %>%
  select(outcome, process,
         b_full = b_interact, p_full = p_interact, pfdr_full = p_interact_fdr,
         nsubj_full = n_subj, nobs_full = n_obs) %>%
  left_join(
    res_sens %>% select(outcome, process,
                        b_sens = b_interact, p_sens = p_interact, pfdr_sens = p_interact_fdr,
                        nsubj_sens = n_subj, nobs_sens = n_obs),
    by = c("outcome","process")
  ) %>%
  mutate(
    db        = b_sens - b_full,
    sign_flip = sign(b_full) != sign(b_sens),
    sig_full  = pfdr_full < .05,
    sig_sens  = pfdr_sens < .05
  )

dir.create("tables", showWarnings = FALSE)
write_csv(cmp, "tables/sensitivity_min30obs.csv")

cat("\n=== Full sample vs >= 30 observations ===\n")
cat(sprintf("Models compared: %d\n", nrow(cmp)))
cat(sprintf("Interactions with FDR p < .05, full: %d | restricted: %d\n",
            sum(cmp$sig_full, na.rm = TRUE), sum(cmp$sig_sens, na.rm = TRUE)))
cat(sprintf("Sign changes of the interaction: %d\n", sum(cmp$sign_flip, na.rm = TRUE)))
cat(sprintf("Correlation of coefficients (full vs restricted): %.3f\n",
            cor(cmp$b_full, cmp$b_sens, use = "complete.obs")))
cat(sprintf("Median absolute change: %.4f | maximum: %.4f\n",
            median(abs(cmp$db), na.rm = TRUE), max(abs(cmp$db), na.rm = TRUE)))
cat("\nRange of the interaction coefficient:\n")
cat(sprintf("  full:       [%.3f, %.3f]\n", min(cmp$b_full, na.rm=T), max(cmp$b_full, na.rm=T)))
cat(sprintf("  restricted: [%.3f, %.3f]\n", min(cmp$b_sens, na.rm=T), max(cmp$b_sens, na.rm=T)))
cat(sprintf("Participants in models, full: %d | restricted: %d\n", max(cmp$nsubj_full), max(cmp$nsubj_sens)))
