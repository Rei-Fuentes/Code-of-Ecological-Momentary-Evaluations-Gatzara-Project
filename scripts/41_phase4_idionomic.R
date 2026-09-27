# =============================================================================
# 41_phase4_idionomic.R · Phase 4: one coefficient per person
# =============================================================================
# What we do here
#   Instead of one average effect, we estimate each participant's own lagged coefficient for every
#   pair: lm(Y(t+1) ~ Y(t) + X(t-1)) on triads of consecutive days, requiring at least 15 triads.
#   We then pool the coefficients of both arms with a random-effects meta-analysis per pair (metafor,
#   REML; FDR across the 24 pairs), and group participants with k-means within each arm (complete
#   24-coefficient profiles, standardized, k from 2 to 5 chosen by mean silhouette, 50 starts).
#   The meta-analysis within each arm, the one in Table 2, is in 42_phase4_meta_by_arm.R.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/idionomic_betas.csv
#         tables/idionomic_meta.csv
#         tables/idionomic_clusters.csv
# Reported in  Results, Phase 4; Tables S5 and S13; Figure S4; Supplementary S8
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(metafor); library(cluster)
})

set.seed(20260519)

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
d  <- df %>% filter(phase == "intervention") %>% arrange(code, date)

outcomes  <- c("global_wellbeing","hedonic_wellbeing","cognitive_wellbeing","negative_emotion_regulation")
processes <- c("emotional_awareness","compassion","self_compassion",
               "positive_emotion_regulation","gratitude","cognitive_reappraisal")

# Participant-specific coefficient of X(t-1), controlling for Y(t).
fit_subject <- function(sub, Y, X) {
  sub <- sub %>%
    arrange(date) %>%
    mutate(
      Y_t   = !!sym(Y),
      Y_tp1 = lead(!!sym(Y), 1),
      X_tm1 = lag(!!sym(X), 1),
      d_prev = as.integer(date - lag(date)),
      d_next = as.integer(lead(date) - date)
    ) %>%
    filter(d_prev == 1, d_next == 1) %>%
    filter(complete.cases(Y_t, Y_tp1, X_tm1))
  if (nrow(sub) < 15) return(c(NA, NA))
  fit <- tryCatch(lm(Y_tp1 ~ Y_t + X_tm1, data = sub), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA))
  s <- summary(fit)$coefficients
  if (!"X_tm1" %in% rownames(s)) return(c(NA, NA))
  c(s["X_tm1","Estimate"], s["X_tm1","Std. Error"])
}

idio_rows <- list(); i <- 1
for (id in unique(d$code)) {
  sub <- d %>% filter(code == id)
  arm <- sub$arm[1]
  for (Y in outcomes) for (X in processes) {
    est <- fit_subject(sub, Y, X)
    idio_rows[[i]] <- tibble(code = id, arm = arm, outcome = Y, process = X,
                             beta = est[1], se = est[2])
    i <- i + 1
  }
}
idio <- bind_rows(idio_rows)
write_csv(idio, "tables/idionomic_betas.csv")

# Random-effects meta-analysis per pair, both arms pooled.
meta_rows <- list(); i <- 1
for (Y in outcomes) for (X in processes) {
  sub <- idio %>% filter(outcome == Y, process == X, !is.na(beta), !is.na(se), se > 0)
  if (nrow(sub) < 5) next
  m <- tryCatch(rma(yi = sub$beta, sei = sub$se, method = "REML"), error = function(e) NULL)
  if (is.null(m)) next
  pi <- predict(m)
  meta_rows[[i]] <- tibble(outcome = Y, process = X, k = m$k,
                           beta_pooled = as.numeric(m$b), se_pooled = m$se,
                           ci_lo = m$ci.lb, ci_hi = m$ci.ub,
                           tau2 = m$tau2, I2 = m$I2,
                           pi_lo = pi$pi.lb, pi_hi = pi$pi.ub,
                           p_value = m$pval)
  i <- i + 1
}
meta <- bind_rows(meta_rows)
if (nrow(meta) > 0) meta <- meta %>% mutate(p_fdr = p.adjust(p_value, method = "BH"))
write_csv(meta, "tables/idionomic_meta.csv")

cat("\n=== Pooled meta-analysis (highest heterogeneity first) ===\n")
print(meta %>% arrange(desc(I2)) %>% head(8) %>%
        select(outcome, process, k, beta_pooled, I2, p_value, p_fdr))

# k-means within arm, k chosen by mean silhouette.
mat <- idio %>%
  mutate(pair = paste(outcome, process, sep = "__")) %>%
  select(code, arm, pair, beta) %>%
  pivot_wider(names_from = pair, values_from = beta) %>%
  filter(if_all(-c(code, arm), ~ !is.na(.)))

cat(sprintf("\nComplete coefficient matrix: %d participants x %d pairs\n",
            nrow(mat), ncol(mat) - 2))

cluster_arm <- function(arm_data, label) {
  X <- as.matrix(arm_data %>% select(-code, -arm))
  X <- scale(X)
  sil <- sapply(2:5, function(k) {
    cl <- kmeans(X, centers = k, nstart = 50)
    mean(silhouette(cl$cluster, dist(X))[, 3])
  })
  best_k <- which.max(sil) + 1
  cl <- kmeans(X, centers = best_k, nstart = 50)
  tibble(code = arm_data$code, arm = label,
         cluster = cl$cluster,
         silhouette_best_k = best_k,
         silhouette_avg = sil[best_k - 1])
}

cluster_assign <- bind_rows(
  cluster_arm(mat %>% filter(arm == "CBWT"), "CBWT"),
  cluster_arm(mat %>% filter(arm == "CBT"), "CBT")
)
write_csv(cluster_assign, "tables/idionomic_clusters.csv")

cat("\n=== Cluster sizes ===\n")
print(cluster_assign %>% group_by(arm, cluster) %>% summarise(n = n(), .groups = "drop") %>%
        ungroup() %>% left_join(
          cluster_assign %>% distinct(arm, silhouette_best_k, silhouette_avg),
          by = "arm"
        ))

cat("\nWritten to tables/.\n")
