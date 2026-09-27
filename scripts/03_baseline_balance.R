# =============================================================================
# 03_baseline_balance.R · Did the arms start from the same place?
# =============================================================================
# What we do here
#   Randomization does not guarantee equal starting points, so we compare the two arms on the
#   baseline-phase prompts of every daily variable. Prompts are nested in people, so we use one
#   multilevel model per variable (value ~ arm + (1 | participant)) instead of a t test on prompts,
#   and we correct the ten p values with the Benjamini-Hochberg false discovery rate.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/baseline_balance.csv
# Reported in  Results, Participant flow; Discussion, Limitations
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(lme4); library(lmerTest) })

# The strength-based arm is the reference, so each coefficient is barrier-based minus strength-based.
d <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE) %>%
  filter(phase == "baseline") %>% mutate(arm = factor(arm, levels = c("CBWT", "CBT")))
V <- c("hedonic_wellbeing","cognitive_wellbeing","global_wellbeing","emotional_awareness","compassion",
       "self_compassion","positive_emotion_regulation","gratitude","negative_emotion_regulation",
       "cognitive_reappraisal")
cat(sprintf("Baseline phase: %d prompts, %d participants\n\n", nrow(d), n_distinct(d$code)))

res <- bind_rows(lapply(V, function(v) {
  m  <- lmer(as.formula(paste(v, "~ arm + (1 | code)")), data = d, REML = TRUE)
  co <- summary(m)$coefficients; row <- grep("arm", rownames(co))
  vc <- as.data.frame(VarCorr(m))
  tibble(variable = v,
         diff_barrier_minus_strength = co[row, "Estimate"],   # barrier-based minus strength-based
         se = co[row, "Std. Error"], p = co[row, "Pr(>|t|)"],
         icc = vc$vcov[1] / sum(vc$vcov))
})) %>% mutate(p_fdr = p.adjust(p, method = "BH"))

write_csv(res, "tables/baseline_balance.csv")
print(res %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf)
items <- res %>% filter(variable != "global_wellbeing")
cat(sprintf("\nItems higher in the barrier-based arm: %d of %d\n",
            sum(items$diff_barrier_minus_strength > 0), nrow(items)))
cat(sprintf("Smallest uncorrected p: %.3f (%s); smallest FDR p: %.3f\n",
            min(res$p), res$variable[which.min(res$p)], min(res$p_fdr)))
