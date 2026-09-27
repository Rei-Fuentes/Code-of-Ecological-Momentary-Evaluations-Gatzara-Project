# =============================================================================
# 12_phase1_sensitivity_specification.R · Phase 1 under two alternative specifications
# =============================================================================
# What we do here
#   We refit the 24 Phase 1 models without the baseline term, in two ways:
#     A. no_baseline_centered: same within-person centering as the main models. This isolates what
#        the baseline adjustment contributes (in these data, almost nothing).
#     B. no_baseline_zwithin: variables standardized within each person (z scores), as in Laicher
#        et al. (2025, Model 2). This was the first specification fitted in the project, and the only
#        one in which two emotional-awareness interactions reach significance.
#
# Input   data/processed/ema_analysis_ready.csv
#         tables/phase1_lagged_mlm.csv
# Output  tables/sensitivity_specification.csv
# Reported in  Method, Analytic Strategy; Results, Phase 1; Discussion
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
d  <- df %>% filter(phase == "intervention")

outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")
vars <- unique(c(outcomes, processes))

cwc <- function(x) x - mean(x, na.rm = TRUE)
zwithin <- function(x) {
  if (sum(!is.na(x)) < 2 || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

fit_all <- function(transform, label) {
  dd <- d %>% arrange(code, date) %>% group_by(code) %>%
    mutate(across(all_of(vars), transform, .names = "{.col}_s")) %>% ungroup()
  res <- list(); i <- 1
  for (Y in outcomes) for (X in processes) {
    dat <- dd %>%
      arrange(code, date) %>% group_by(code) %>%
      mutate(
        Y_t    = !!sym(paste0(Y, "_s")),
        Y_tp1  = lead(!!sym(paste0(Y, "_s")), 1),
        X_t    = !!sym(paste0(X, "_s")),
        X_tm1  = lag(!!sym(paste0(X, "_s")), 1),
        d_prev = as.integer(date - lag(date)),
        d_next = as.integer(lead(date) - date)
      ) %>% ungroup() %>%
      filter(d_prev == 1, d_next == 1) %>%
      filter(!is.na(Y_t), !is.na(Y_tp1), !is.na(X_t), !is.na(X_tm1)) %>%
      mutate(arm = factor(arm, levels = c("CBT","CBWT")))
    if (nrow(dat) < 50) next
    fit <- tryCatch(lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + (1 | code), data = dat, REML = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)$coefficients
    b <- if ("X_tm1:armCBWT" %in% rownames(s)) s["X_tm1:armCBWT", c("Estimate","Std. Error","Pr(>|t|)")] else c(NA,NA,NA)
    res[[i]] <- tibble(specification = label, outcome = Y, process = X,
                       n_obs = nrow(dat), n_subj = n_distinct(dat$code),
                       b_interact = b[1], se_interact = b[2], p_interact = b[3])
    i <- i + 1
  }
  bind_rows(res) %>% mutate(p_interact_fdr = p.adjust(p_interact, method = "BH"))
}

out <- bind_rows(fit_all(cwc, "no_baseline_centered"),
                 fit_all(zwithin, "no_baseline_zwithin")) %>%
  arrange(specification, outcome, process)
write_csv(out, "tables/sensitivity_specification.csv")

primary <- read_csv("tables/phase1_lagged_mlm.csv", show_col_types = FALSE) %>%
  select(outcome, process, b_primary = b_interact)
for (sp in unique(out$specification)) {
  o <- out %>% filter(specification == sp) %>% left_join(primary, by = c("outcome","process"))
  cat(sprintf("\n=== %s ===\nInteractions with FDR p < .05: %d | max |b - b_primary| = %.4f\n",
              sp, sum(o$p_interact_fdr < .05), max(abs(o$b_interact - o$b_primary))))
  print(o %>% arrange(p_interact_fdr) %>% head(3) %>%
          select(outcome, process, b_interact, p_interact, p_interact_fdr))
}
