# =============================================================================
# 18_phase1_power_simr.R · How large an interaction could we detect?
# =============================================================================
# What we do here
#   A null result is only informative if we know what the design could detect. We take one
#   representative Phase 1 model (reappraisal predicting negative-emotion regulation, the pair with
#   the largest interaction), set its interaction to a range of values, and simulate new data 1,000
#   times per value with simr. Power is the share of simulations with p < .05. Coefficients are in
#   scale points, because the daily variables are centered within person on their 1-7 scale.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/sensitivity_power_simr.csv
# Reported in  Method, Sample size and power; Supplementary S4
#
# About 30 minutes.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(lme4); library(lmerTest); library(simr)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)

# Representative pair
Y <- "negative_emotion_regulation"
X <- "cognitive_reappraisal"
vars <- c(Y, X)

baseline <- df %>% filter(phase == "baseline") %>%
  group_by(code) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_base"),
            .groups = "drop")
cwc <- function(x) x - mean(x, na.rm = TRUE)

d <- df %>% filter(phase == "intervention") %>%
  left_join(baseline, by = "code") %>%
  arrange(code, date) %>% group_by(code) %>%
  mutate(across(all_of(vars), cwc, .names = "{.col}_c")) %>% ungroup()

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

cat(sprintf("Data: %d triads, %d participants\n",
            nrow(dat), n_distinct(dat$code)))

fit <- lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + Y_base + (1 | code),
            data = dat, REML = FALSE)
cat("\nFitted model:\n")
print(round(summary(fit)$coefficients, 4))

# Interaction values to simulate
target_effects <- c(0.03, 0.05, 0.07, 0.10, 0.13)
n_sims <- 1000  # simulations per effect size

cat(sprintf("\nSimulating power for X_tm1:armCBWT, %d simulations per value...\n", n_sims))
results <- data.frame(effect = target_effects, power = NA_real_,
                      ci_lo = NA_real_, ci_hi = NA_real_)

for (i in seq_along(target_effects)) {
  e <- target_effects[i]
  fit_sim <- fit
  fixef(fit_sim)["X_tm1:armCBWT"] <- e
  ps <- powerSim(fit_sim, test = fixed("X_tm1:armCBWT", method = "t"),
                 nsim = n_sims, progress = FALSE, alpha = .05)
  s <- summary(ps)
  results$power[i] <- s$mean
  results$ci_lo[i] <- s$lower
  results$ci_hi[i] <- s$upper
  cat(sprintf("  effect=%.3f  -> power=%.2f%% [95%% CI %.1f, %.1f]\n",
              e, 100*s$mean, 100*s$lower, 100*s$upper))
}

write_csv(results, "tables/sensitivity_power_simr.csv")
cat("\nWritten: tables/sensitivity_power_simr.csv\n")

# Smallest simulated value with >= 80% power
above <- which(results$power >= 0.80)
if (length(above) > 0) {
  min_detectable <- results$effect[min(above)]
  cat(sprintf("\nSmallest simulated value with >= 80%% power: |b| >= %.3f\n",
              min_detectable))
} else {
  cat("\nNo simulated value reaches 80% power.\n")
}
