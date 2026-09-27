# =============================================================================
# 33_phase3_minsize_sensitivity.R · Phase 3: does min.size change the picture?
# =============================================================================
# What we do here
#   We repeat the change-point analysis with min.size of 10, 14 and 21.
#     A. As first run during the project's statistical audit: all change-points (gains or not) in six
#        key variables, with 200 permutations. Kept for the record.
#     B. The full gain rule of 31_phase3_changepoints.R on the ten variables, with 1,000
#        permutations. We reset the seed before each min.size, so min.size = 14 reproduces the main
#        analysis exactly. This is the version reported.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/changepoint_minsize_sensitivity.csv
# Reported in  Supplementary S7
#
# About 15 minutes.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(ecp) })

d <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE) %>%
  filter(phase == "intervention") %>% arrange(code, date)

# ---- A. Change-point counts, six key variables, R = 200 ------------------------------------
vars_a <- c("global_wellbeing","self_compassion","emotional_awareness",
            "compassion","negative_emotion_regulation","cognitive_reappraisal")
count_cps <- function(min_size) {
  total <- 0
  for (id in unique(d$code)) {
    sub <- d %>% filter(code == id)
    if (nrow(sub) < min_size * 2) next
    for (v in vars_a) {
      y <- sub[[v]]
      if (sum(!is.na(y)) < min_size * 2) next
      fit <- tryCatch(e.divisive(matrix(y, ncol = 1), sig.lvl = .005, R = 200,
                                 min.size = min_size, alpha = 2), error = function(e) NULL)
      if (is.null(fit)) next
      cps <- fit$estimates; cps <- cps[cps > 1 & cps < length(y)]
      total <- total + length(cps)
    }
  }
  total
}

# ---- B. Gains under the full rule, ten variables, R = 1,000 --------------------------------
vars_b <- c("hedonic_wellbeing","cognitive_wellbeing","global_wellbeing","emotional_awareness",
            "compassion","self_compassion","positive_emotion_regulation",
            "gratitude","negative_emotion_regulation","cognitive_reappraisal")
trajectory_improves <- function(y) {
  if (length(y) < 14) return(FALSE)
  fit <- tryCatch(
    e.divisive(matrix(y, ncol = 1), sig.lvl = .05, R = 199, k = 1,
               min.size = max(7, floor(length(y)/3)), alpha = 2),
    error = function(e) NULL
  )
  if (is.null(fit) || length(fit$estimates) < 3) return(FALSE)
  cp <- fit$estimates[2]
  if (is.na(cp) || cp < 2 || cp >= length(y)) return(FALSE)
  mean(y[cp:length(y)], na.rm = TRUE) > mean(y[1:(cp-1)], na.rm = TRUE)
}
count_gains <- function(min_size) {
  n_cp <- 0; n_gain <- 0
  for (id in unique(d$code)) {
    sub <- d %>% filter(code == id)
    for (v in vars_b) {
      y <- sub[[v]]
      if (sum(!is.na(y)) < 28) next    # same series as the main analysis
      if (length(y) < min_size * 2) next
      fit <- tryCatch(e.divisive(matrix(y, ncol = 1), sig.lvl = .005, R = 1000, k = NULL,
                                 min.size = min_size, alpha = 2), error = function(e) NULL)
      if (is.null(fit)) next
      cps <- fit$estimates; cps <- cps[cps > 1 & cps < length(y)]
      if (length(cps) == 0) next
      up <- trajectory_improves(y)
      for (cp in cps) {
        w_before <- y[max(1, cp - 14):(cp - 1)]; w_after <- y[cp:min(length(y), cp + 13)]
        delta <- mean(w_after, na.rm = TRUE) - mean(w_before, na.rm = TRUE)
        n_cp <- n_cp + 1; n_gain <- n_gain + (!is.na(delta) && delta > 0 && up)
      }
    }
  }
  c(n_cp, n_gain)
}

out <- list()
for (ms in c(14, 10, 21)) {
  set.seed(20260519); a <- count_cps(ms)
  set.seed(20260519); b <- count_gains(ms)
  out[[length(out) + 1]] <- tibble(min_size = ms, A_changepoints_6vars_R200 = a,
                                   B_changepoints_10vars = b[1], B_gains_10vars = b[2])
  cat(sprintf("min.size = %d: A %d change-points | B %d change-points, %d gains\n", ms, a, b[1], b[2]))
}
res <- bind_rows(out) %>% arrange(min_size)
write_csv(res, "tables/changepoint_minsize_sensitivity.csv")
print(res)
