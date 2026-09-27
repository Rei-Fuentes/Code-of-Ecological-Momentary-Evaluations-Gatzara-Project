# =============================================================================
# 17_phase1_slopes_by_arm.R · Phase 1 slopes within each arm
# =============================================================================
# What we do here
#   For the figure we need the lagged slope in each arm, not only their difference. The barrier-based
#   slope is b(X(t-1)); the strength-based slope is b(X(t-1)) + b(X(t-1):armCBWT), and its standard
#   error must include the covariance of the two terms: sqrt(V11 + V22 + 2 * V12).
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/lagged_slopes_by_arm.csv
# Reported in  Figure S2
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(lme4); library(lmerTest) })
df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion","positive_emotion_regulation","gratitude","cognitive_reappraisal")
vars <- unique(c(outcomes, processes))
baseline <- df %>% filter(phase == "baseline") %>% group_by(code) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_base"), .groups = "drop")
d <- df %>% filter(phase == "intervention") %>% left_join(baseline, by = "code")
cwc <- function(x) x - mean(x, na.rm = TRUE)
d <- d %>% arrange(code, date) %>% group_by(code) %>% mutate(across(all_of(vars), cwc, .names = "{.col}_c")) %>% ungroup()
out <- list()
for (Y in outcomes) for (X in processes) {
  dat <- d %>% arrange(code, date) %>% group_by(code) %>%
    mutate(Y_t = !!sym(paste0(Y,"_c")), Y_tp1 = lead(!!sym(paste0(Y,"_c")), 1),
           X_t = !!sym(paste0(X,"_c")), X_tm1 = lag(!!sym(paste0(X,"_c")), 1), Y_base = !!sym(paste0(Y,"_base")),
           d_prev = as.integer(date - lag(date)), d_next = as.integer(lead(date) - date)) %>% ungroup() %>%
    filter(d_prev == 1, d_next == 1) %>% filter(!is.na(Y_t), !is.na(Y_tp1), !is.na(X_t), !is.na(X_tm1), !is.na(Y_base)) %>%
    mutate(arm = factor(arm, levels = c("CBT","CBWT")))
  f <- lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + Y_base + (1 | code), data = dat, REML = TRUE)
  b <- fixef(f); V <- as.matrix(vcov(f)); i1 <- "X_tm1"; i2 <- "X_tm1:armCBWT"
  out[[length(out)+1]] <- tibble(outcome = Y, process = X,
    b_barrier = b[i1], se_barrier = sqrt(V[i1,i1]),
    b_strength = b[i1] + b[i2], se_strength = sqrt(V[i1,i1] + V[i2,i2] + 2*V[i1,i2]),
    se_strength_no_cov = sqrt(V[i1,i1] + V[i2,i2]))
}
res <- bind_rows(out); write_csv(res, "tables/lagged_slopes_by_arm.csv")
r <- res$se_strength / res$se_strength_no_cov
cat(sprintf("SE with / without the covariance term: median %.3f, range [%.3f, %.3f]\n", median(r), min(r), max(r)))
ci0 <- with(res, (b_strength - 1.96*se_strength_no_cov) > 0 | (b_strength + 1.96*se_strength_no_cov) < 0)
ci1 <- with(res, (b_strength - 1.96*se_strength) > 0 | (b_strength + 1.96*se_strength) < 0)
cat(sprintf("Strength-based 95%% CIs excluding 0: without covariance %d, with covariance %d (of 24)\n", sum(ci0), sum(ci1)))
