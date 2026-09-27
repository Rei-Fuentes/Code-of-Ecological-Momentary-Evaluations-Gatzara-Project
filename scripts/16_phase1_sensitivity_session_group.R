# =============================================================================
# 16_phase1_sensitivity_session_group.R · Does the session group matter?
# =============================================================================
# What we do here
#   Participants chose one of three weekdays for their group, which gives six session groups nested
#   in the two arms. The main models ignore this grouping, so we check it in two ways:
#     A. How much of the variance in daily scores is due to session group and to participant:
#        value ~ arm + (1 | session_group) + (1 | participant).
#     B. The 24 Phase 1 models refitted with a random intercept for session group.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/sensitivity_group_icc.csv
#         tables/sensitivity_group_phase1.csv
# Reported in  Method, Analytic Strategy; Results, Phase 1; Table S14
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest)
})
set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")
vars <- unique(c(outcomes, processes))

# ---- A. Session-group ICC, net of arm ------------------------------------
dur <- df %>% filter(phase == "intervention") %>%
  mutate(arm = factor(arm, levels = c("CBT","CBWT")))
icc_rows <- lapply(vars, function(v) {
  f <- lmer(as.formula(paste(v, "~ arm + (1 | session_group) + (1 | code)")), data = dur, REML = TRUE)
  vc <- as.data.frame(VarCorr(f))
  g <- vc$vcov[vc$grp == "session_group"]; p <- vc$vcov[vc$grp == "code"]; r <- vc$vcov[vc$grp == "Residual"]
  tibble(variable = v, var_group = g, var_participant = p, var_residual = r,
         icc_group = g / (g + p + r), icc_participant = p / (g + p + r),
         singular = isSingular(f))
})
icc <- bind_rows(icc_rows)
write_csv(icc, "tables/sensitivity_group_icc.csv")

# ---- B. Phase 1 with (1 | session_group) -----------------------------------------
baseline <- df %>% filter(phase == "baseline") %>% group_by(code) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_base"), .groups = "drop")
d <- df %>% filter(phase == "intervention") %>% left_join(baseline, by = "code")
cwc <- function(x) x - mean(x, na.rm = TRUE)
d <- d %>% arrange(code, date) %>% group_by(code) %>%
  mutate(across(all_of(vars), cwc, .names = "{.col}_c")) %>% ungroup()

res <- list(); i <- 1
for (Y in outcomes) for (X in processes) {
  dat <- d %>% arrange(code, date) %>% group_by(code) %>%
    mutate(Y_t = !!sym(paste0(Y,"_c")), Y_tp1 = lead(!!sym(paste0(Y,"_c")), 1),
           X_t = !!sym(paste0(X,"_c")), X_tm1 = lag(!!sym(paste0(X,"_c")), 1),
           Y_base = !!sym(paste0(Y,"_base")),
           d_prev = as.integer(date - lag(date)),
           d_next = as.integer(lead(date) - date)) %>% ungroup() %>%
    filter(d_prev == 1, d_next == 1) %>%
    filter(!is.na(Y_t), !is.na(Y_tp1), !is.na(X_t), !is.na(X_tm1), !is.na(Y_base)) %>%
    mutate(arm = factor(arm, levels = c("CBT","CBWT")))
  f0 <- lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + Y_base + (1 | code), data = dat, REML = TRUE)
  f1 <- lmer(Y_tp1 ~ Y_t + X_t + X_tm1 * arm + Y_base + (1 | session_group) + (1 | code), data = dat, REML = TRUE)
  s0 <- summary(f0)$coefficients; s1 <- summary(f1)$coefficients
  vc1 <- as.data.frame(VarCorr(f1))
  res[[i]] <- tibble(outcome = Y, process = X, n_subj = n_distinct(dat$code),
    b_orig = s0["X_tm1:armCBWT","Estimate"], p_orig = s0["X_tm1:armCBWT","Pr(>|t|)"],
    b_group = s1["X_tm1:armCBWT","Estimate"], p_group = s1["X_tm1:armCBWT","Pr(>|t|)"],
    var_group = vc1$vcov[vc1$grp == "session_group"], singular_group = isSingular(f1))
  i <- i + 1
}
cmp <- bind_rows(res) %>%
  mutate(pfdr_orig = p.adjust(p_orig, "BH"), pfdr_group = p.adjust(p_group, "BH"),
         db = b_group - b_orig, sign_flip = sign(b_group) != sign(b_orig))
write_csv(cmp, "tables/sensitivity_group_phase1.csv")

cat("=== A. Session-group ICC (net of arm), intervention phase ===\n")
print(icc %>% mutate(across(where(is.numeric), ~ round(., 4))) %>%
        select(variable, icc_group, icc_participant, singular))
cat(sprintf("\nSession-group ICC: median %.4f, maximum %.4f\n", median(icc$icc_group), max(icc$icc_group)))
cat(sprintf("Participant ICC: median %.3f\n", median(icc$icc_participant)))

cat("\n=== B. Phase 1 with (1 | session_group) vs original ===\n")
cat(sprintf("Interactions with FDR p < .05, original: %d, with group: %d\n",
            sum(cmp$pfdr_orig < .05), sum(cmp$pfdr_group < .05)))
cat(sprintf("Sign changes: %d\n", sum(cmp$sign_flip)))
cat(sprintf("Correlation of coefficients: %.4f\n", cor(cmp$b_orig, cmp$b_group)))
cat(sprintf("|change| median %.5f, maximum %.5f\n", median(abs(cmp$db)), max(abs(cmp$db))))
cat(sprintf("Models with session-group variance = 0 (singular): %d of %d\n",
            sum(cmp$singular_group), nrow(cmp)))
cat("\nLowest-p interaction in both specifications:\n")
print(cmp %>% arrange(p_orig) %>% head(3) %>%
        select(outcome, process, b_orig, pfdr_orig, b_group, pfdr_group) %>%
        mutate(across(where(is.numeric), ~ round(., 4))))
