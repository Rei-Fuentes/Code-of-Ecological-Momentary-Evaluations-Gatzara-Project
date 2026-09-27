# =============================================================================
# 32_phase3_gain_rate_contrast.R · Phase 3: are gains more frequent in one arm?
# =============================================================================
# What we do here
#   We count gains per analysable participant (at least 28 intervention-phase prompts), including
#   those with none, and compare arms with the rate ratio. Its 95% interval comes from 10,000
#   bootstrap resamples of participants within arm; its p value from 10,000 permutations of the arm
#   labels (we report the two-sided value and also save the one-sided one). We also compare the share
#   of participants with at least one gain (chi-square test with continuity correction).
#
# Input   data/processed/ema_analysis_ready.csv
#         tables/changepoints_all.csv
# Output  tables/changepoint_rate_contrast.csv
# Reported in  Abstract; Results, Phase 3
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })

set.seed(20260519)
B <- 10000

df  <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
cps <- read_csv("tables/changepoints_all.csv", show_col_types = FALSE)

analysable <- df %>% filter(phase == "intervention") %>% count(code, arm, name = "n_obs") %>%
  filter(n_obs >= 28) %>% select(code, arm)

gains <- cps %>% filter(is_gain)
per_person <- analysable %>%
  left_join(gains %>% count(code, name = "n_gains"), by = "code") %>%
  mutate(n_gains = coalesce(n_gains, 0L))

g_s <- per_person$n_gains[per_person$arm == "CBWT"]
g_b <- per_person$n_gains[per_person$arm == "CBT"]
rr  <- function(a, b) mean(a) / mean(b)
rr_obs <- rr(g_s, g_b)

boot <- replicate(B, rr(sample(g_s, replace = TRUE), sample(g_b, replace = TRUE)))
all_g <- per_person$n_gains; n_s <- length(g_s)
perm <- replicate(B, { i <- sample(length(all_g), n_s); rr(all_g[i], all_g[-i]) })

any_gain <- table(factor(per_person$arm, levels = c("CBWT","CBT")), per_person$n_gains > 0)
chi <- suppressWarnings(chisq.test(any_gain, correct = TRUE))

pairs_all <- gains %>% distinct(code, variable)
res <- tibble(
  analysable_strength = length(g_s), analysable_barrier = length(g_b),
  gains_total = nrow(gains), gains_strength = sum(g_s), gains_barrier = sum(g_b),
  participant_variable_pairs_with_gain = nrow(pairs_all),
  pairs_with_gain_strength = nrow(pairs_all %>% semi_join(analysable %>% filter(arm == "CBWT"), by = "code")),
  pairs_with_gain_barrier  = nrow(pairs_all %>% semi_join(analysable %>% filter(arm == "CBT"), by = "code")),
  rate_strength = mean(g_s), rate_barrier = mean(g_b), rate_ratio = rr_obs,
  boot_ci_lo = unname(quantile(boot, .025)), boot_ci_hi = unname(quantile(boot, .975)),
  perm_p_one_sided = mean(perm >= rr_obs),
  perm_p_two_sided = mean(abs(log(perm)) >= abs(log(rr_obs))),
  zero_gain_strength = sum(g_s == 0), zero_gain_barrier = sum(g_b == 0),
  pct_any_gain_strength = 100 * mean(g_s > 0), pct_any_gain_barrier = 100 * mean(g_b > 0),
  chisq_p_any_gain = chi$p.value, B = B)
write_csv(res, "tables/changepoint_rate_contrast.csv")
print(t(res %>% mutate(across(where(is.numeric), ~ round(., 3)))))
