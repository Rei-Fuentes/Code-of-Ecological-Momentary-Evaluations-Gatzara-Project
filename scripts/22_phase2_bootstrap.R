# =============================================================================
# 22_phase2_bootstrap.R · Phase 2: how stable are the networks?
# =============================================================================
# What we do here
#   We resample participants with replacement within each arm (1,000 times per arm), refit the
#   network each time and keep the strength of every node and the global strength of each network.
#   That gives percentile intervals and tells us how often each node ranks first or in the top three.
#   We draw all resamples first and fit them afterwards. mlVAR does not use random numbers, so the
#   fits can run on several cores (N_CORES) and the result is the same as fitting them one by one.
#
# Input   data/processed/ema_analysis_ready.csv
# Output  tables/mlvar_centrality_bootstrap_B1000_raw.csv
#         tables/mlvar_centrality_bootstrap_B1000_summary.csv
#         tables/mlvar_global_strength_B1000.csv
# Reported in  Tables S2, S3 and S6; Supplementary S6
#
# Long run: about 7 hours on one core, about 1 hour with N_CORES=8.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(mlVAR); library(qgraph)
})

set.seed(20260519)
B <- 1000

df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)
d  <- df %>% filter(phase == "intervention")

vars <- c("hedonic_wellbeing","cognitive_wellbeing","emotional_awareness",
          "compassion","self_compassion","positive_emotion_regulation",
          "gratitude","negative_emotion_regulation","cognitive_reappraisal")

prep_arm <- function(arm_code) {
  d %>% filter(arm == arm_code) %>%
    group_by(code) %>% filter(n() >= 20) %>% ungroup()
}

run_one <- function(arm_data) {
  fit <- mlVAR(as.data.frame(arm_data), vars = vars, idvar = "code",
               lags = 1, estimator = "lmer",
               temporal = "orthogonal", contemporaneous = "orthogonal",
               verbose = FALSE)
  ct <- centrality(getNet(fit, "temporal", nonsig = "show"))
  cc <- centrality(getNet(fit, "contemporaneous", rule = "and", nonsig = "show"))
  gs_t <- sum(abs(getNet(fit, "temporal", nonsig = "show")))
  gs_c <- sum(abs(getNet(fit, "contemporaneous", rule = "and", nonsig = "show")))
  bind_rows(
    tibble(node = vars, network = "temporal",
           in_strength = ct$InDegree, out_strength = ct$OutDegree,
           strength = ct$InDegree + ct$OutDegree, global_strength = gs_t),
    tibble(node = vars, network = "contemporaneous",
           in_strength = cc$InDegree, out_strength = cc$OutDegree,
           strength = cc$InDegree + cc$OutDegree, global_strength = gs_c)
  )
}

# All resamples are drawn first, in the order CBWT then CBT. mlVAR does not use the random
# number generator, so fitting the resamples in parallel (N_CORES > 1, via fork) gives the
# same results as fitting them one after the other.
N_CORES <- as.integer(Sys.getenv("N_CORES", "1"))

draw_resamples <- function(arm_code) {
  ids <- unique(prep_arm(arm_code)$code)
  lapply(seq_len(B), function(b) sample(ids, length(ids), replace = TRUE))
}
draws <- list(CBWT = draw_resamples("CBWT"), CBT = draw_resamples("CBT"))

bootstrap_arm <- function(arm_code) {
  base <- prep_arm(arm_code)
  cat(sprintf("[%s] Bootstrap %s - N=%d, B=%d, cores=%d\n", format(Sys.time()), arm_code,
              n_distinct(base$code), B, N_CORES))
  one <- function(b) {
    idx <- draws[[arm_code]][[b]]
    boot <- bind_rows(lapply(seq_along(idx), function(i) {
      sub <- base %>% filter(code == idx[i]); sub$code <- paste0(idx[i], "_b", i); sub
    }))
    res <- tryCatch(run_one(boot), error = function(e) NULL)
    if (!is.null(res)) res$boot <- b
    res
  }
  out <- if (N_CORES > 1) parallel::mclapply(seq_len(B), one, mc.cores = N_CORES) else lapply(seq_len(B), one)
  ok <- !vapply(out, is.null, logical(1))
  cat(sprintf("[%s] %s done: %d/%d valid iterations\n", format(Sys.time()), arm_code, sum(ok), B))
  bind_rows(out[ok]) %>% mutate(arm = arm_code)
}

raw <- bind_rows(bootstrap_arm("CBWT"), bootstrap_arm("CBT"))
write_csv(raw, "tables/mlvar_centrality_bootstrap_B1000_raw.csv")

# Per node: bootstrap mean and 95% percentile interval of in-, out- and total strength.
summ <- raw %>%
  pivot_longer(c(in_strength, out_strength, strength), names_to = "index", values_to = "v") %>%
  group_by(arm, network, node, index) %>%
  summarise(mean = mean(v), lo = quantile(v, .025), hi = quantile(v, .975),
            n_boot = n(), .groups = "drop")
write_csv(summ, "tables/mlvar_centrality_bootstrap_B1000_summary.csv")

gs <- raw %>% distinct(arm, network, boot, global_strength) %>%
  group_by(arm, network) %>%
  summarise(mean = mean(global_strength), lo = quantile(global_strength, .025),
            hi = quantile(global_strength, .975), n_boot = n(), .groups = "drop")
write_csv(gs, "tables/mlvar_global_strength_B1000.csv")

cat("\nWritten to tables/.\n")
print(summ %>% filter(network == "temporal", index == "out_strength") %>% arrange(arm, desc(mean)))
