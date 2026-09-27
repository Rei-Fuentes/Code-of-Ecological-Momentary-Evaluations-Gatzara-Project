# =============================================================================
# 11_phase1_lagged_mlm.R · Phase 1: do the daily links differ between arms?
# =============================================================================
# What we do here
#   This is the between-arm test of the paper. For each of the 24 process-outcome pairs we ask
#   whether yesterday's process predicts tomorrow's outcome more strongly in one arm than in the
#   other, using triads of three consecutive days in the intervention phase:
#     Y(t+1) ~ Y(t) + X(t) + X(t-1) * arm + baseline + (1 | participant)
#   Variables are centered within each person; baseline is the person's mean of the outcome in the
#   baseline phase. The barrier-based arm (CBT) is the reference, so X(t-1):armCBWT is the
#   difference in slope. The 24 interaction tests are corrected with the false discovery rate.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/phase1_lagged_mlm.csv
# Reported in  Method, Analytic Strategy; Results, Phase 1; Table S1
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(lme4); library(lmerTest); library(MuMIn)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)

outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")
vars <- unique(c(outcomes, processes))

# Baseline level per participant (mean in the baseline phase).
baseline <- df %>%
  filter(phase == "baseline") %>%
  group_by(code) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_base"),
            .groups = "drop")

d <- df %>% filter(phase == "intervention") %>%
  left_join(baseline, by = "code")

# Within-person centering.
cwc <- function(x) x - mean(x, na.rm = TRUE)
d <- d %>% arrange(code, date) %>% group_by(code) %>%
  mutate(across(all_of(vars), cwc, .names = "{.col}_c")) %>%
  ungroup()

results <- list(); i <- 1
for (Y in outcomes) for (X in processes) {
  dat <- d %>%
    arrange(code, date) %>%
    group_by(code) %>%
    mutate(
      Y_t   = !!sym(paste0(Y,"_c")),
      Y_tp1 = lead(!!sym(paste0(Y,"_c")), 1),
      X_t   = !!sym(paste0(X,"_c")),
      X_tm1 = lag(!!sym(paste0(X,"_c")), 1),
      Y_base = !!sym(paste0(Y,"_base")),
      d_prev = as.integer(date - lag(date)),
      d_next = as.integer(lead(date) - date)
    ) %>%
    ungroup() %>%
    filter(d_prev == 1, d_next == 1) %>%
    filter(!is.na(Y_t), !is.na(Y_tp1), !is.na(X_t), !is.na(X_tm1), !is.na(Y_base)) %>%
    mutate(arm = factor(arm, levels = c("CBT","CBWT")))
  if (nrow(dat) < 50) next

  fit <- tryCatch(
    lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + Y_base + (1 | code), data = dat, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) next

  s   <- summary(fit)$coefficients
  vc  <- as.data.frame(VarCorr(fit))
  r2  <- tryCatch(r.squaredGLMM(fit), error = function(e) c(R2m=NA, R2c=NA))
  pull <- function(row) if (row %in% rownames(s)) s[row, c("Estimate","Std. Error","Pr(>|t|)")] else c(NA,NA,NA)

  beta_xtm1     <- pull("X_tm1")
  beta_interact <- pull("X_tm1:armCBWT")

  results[[i]] <- tibble(
    outcome = Y, process = X, n_obs = nrow(dat), n_subj = n_distinct(dat$code),
    b_X_tm1    = beta_xtm1[1],    se_X_tm1    = beta_xtm1[2],    p_X_tm1    = beta_xtm1[3],
    b_interact = beta_interact[1], se_interact = beta_interact[2], p_interact = beta_interact[3],
    var_intercept = vc$vcov[vc$grp == "code"],
    var_residual  = vc$vcov[vc$grp == "Residual"],
    icc = vc$vcov[vc$grp == "code"] / (vc$vcov[vc$grp == "code"] + vc$vcov[vc$grp == "Residual"]),
    singular = isSingular(fit),
    R2_marg = r2[1], R2_cond = r2[2]
  )
  i <- i + 1
}

res <- bind_rows(results) %>%
  mutate(
    p_X_tm1_fdr    = p.adjust(p_X_tm1,    method = "BH"),
    p_interact_fdr = p.adjust(p_interact, method = "BH"),
    sig_interact_fdr = p_interact_fdr < .05
  ) %>%
  arrange(outcome, process)

write_csv(res, "tables/phase1_lagged_mlm.csv")

cat("\n=== Phase 1 lagged multilevel models ===\n")
cat(sprintf("Models: %d  |  Singular fits: %d  |  Median ICC: %.3f\n",
            nrow(res), sum(res$singular), median(res$icc, na.rm=TRUE)))
cat(sprintf("Interactions with FDR p < .05: %d\n", sum(res$sig_interact_fdr, na.rm=TRUE)))
print(res %>% filter(sig_interact_fdr) %>%
        select(outcome, process, b_X_tm1, b_interact, p_interact, p_interact_fdr))

cat("\nInteractions ordered by p:\n")
print(res %>% arrange(p_interact) %>% head(8) %>%
        select(outcome, process, b_interact, p_interact, p_interact_fdr))

cat(sprintf("\nParticipants contributing triads: %s | triads per model: %s\n",
            paste(unique(res$n_subj), collapse = "/"), paste(unique(res$n_obs), collapse = "/")))
