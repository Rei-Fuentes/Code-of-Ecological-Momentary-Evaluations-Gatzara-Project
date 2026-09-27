# =============================================================================
# 42_phase4_meta_by_arm.R · Phase 4: pooled person-level effects within each arm
# =============================================================================
# What we do here
#   We pool the participant-specific coefficients of 41_phase4_idionomic.R separately in each arm
#   (random-effects meta-analysis, REML) and correct the 24 pairs with the FDR within arm. A pair that
#   is significant in one arm and not in the other does not show that the arms differ (Gelman & Stern,
#   2006); that formal test is the Phase 1 interaction.
#
# Input   tables/idionomic_betas.csv
# Output  tables/idionomic_meta_by_arm.csv
# Reported in  Abstract; Results, Phase 4; Table 2
# =============================================================================

suppressPackageStartupMessages({library(dplyr);library(readr);library(metafor)})
idio <- read_csv("tables/idionomic_betas.csv", show_col_types=FALSE)
outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion","positive_emotion_regulation","gratitude","cognitive_reappraisal")
res <- list(); i <- 1
for (br in c("CBWT","CBT")) for (Y in outcomes) for (X in processes) {
  sub <- idio %>% filter(arm==br, outcome==Y, process==X, !is.na(beta), !is.na(se), se>0)
  if (nrow(sub) < 5) next
  m <- tryCatch(rma(yi=sub$beta, sei=sub$se, method="REML"), error=function(e) NULL)
  if (is.null(m)) next
  res[[i]] <- tibble(arm=br, outcome=Y, process=X, k=m$k, beta=as.numeric(m$b),
                     ci_lo=m$ci.lb, ci_hi=m$ci.ub, I2=m$I2, p=m$pval); i <- i+1
}
out <- bind_rows(res) %>% group_by(arm) %>% mutate(p_fdr=p.adjust(p, method="BH")) %>% ungroup()
write_csv(out, "tables/idionomic_meta_by_arm.csv")
cat(sprintf("\nPairs with FDR p < .05: CBWT %d/%d | CBT %d/%d\n",
  sum(out$p_fdr<.05 & out$arm=="CBWT"), sum(out$arm=="CBWT"),
  sum(out$p_fdr<.05 & out$arm=="CBT"), sum(out$arm=="CBT")))
cat("\n--- Selected processes by arm ---\n")
out %>% filter(process %in% c("compassion","self_compassion","cognitive_reappraisal","gratitude"),
               outcome %in% c("global_wellbeing","cognitive_wellbeing")) %>%
  arrange(process, outcome, arm) %>%
  mutate(across(c(beta,ci_lo,ci_hi,p_fdr), ~round(.,3))) %>%
  select(process,outcome,arm,k,beta,ci_lo,ci_hi,p_fdr) %>% print(n=40)

cat(sprintf("\nMedian I2: CBWT %.1f%% | CBT %.1f%%\n",
            median(out$I2[out$arm == "CBWT"]), median(out$I2[out$arm == "CBT"])))
cat(sprintf("k range: CBWT %d-%d | CBT %d-%d\n",
            min(out$k[out$arm == "CBWT"]), max(out$k[out$arm == "CBWT"]),
            min(out$k[out$arm == "CBT"]), max(out$k[out$arm == "CBT"])))
