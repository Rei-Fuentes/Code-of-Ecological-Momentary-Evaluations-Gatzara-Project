# =============================================================================
# 51_figures.R · Data figures
# =============================================================================
# What we do here
#   We draw the figures that come from data:
#     Figure 2   temporal networks by arm             (21_phase2_mlvar.R)
#     Figure S1  weekly means of the daily items      (01_preprocess.R)
#     Figure S2  lagged slopes within each arm        (17_phase1_slopes_by_arm.R)
#     Figure S3  change-points of one participant     (31_phase3_changepoints.R)
#     Figure S4  person-level coefficients by cluster (41_phase4_idionomic.R)
#   Figure 1 (participant flow) is drawn by 52_figure1_consort.py. Figure 3 is a conceptual diagram
#   and does not come from data.
#
# Input   outputs of the scripts above
# Output  figures/*.png
# Reported in  Figures 2 and S1-S4
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(qgraph); library(forcats); library(stringr); library(mlVAR)
})

set.seed(20260519)

FONT_FAMILY <- "Helvetica"

dir.create("figures", showWarnings = FALSE)
df <- read_csv("data/processed/ema_analysis_ready.csv", show_col_types = FALSE)

ARM_STR <- "Strength-based"; ARM_BAR <- "Barrier-based"
pretty <- function(x) x %>% str_replace_all("_", " ") %>%
  str_replace("wellbeing", "well-being") %>% str_replace("self compassion", "self-compassion")
d  <- df %>% filter(phase == "intervention") %>%
  mutate(arm = ifelse(arm == "CBT", ARM_BAR, ARM_STR))

vars <- c("hedonic_wellbeing","cognitive_wellbeing","emotional_awareness",
          "compassion","self_compassion","positive_emotion_regulation",
          "gratitude","negative_emotion_regulation","cognitive_reappraisal")

# Node abbreviations for the network figure.
abbr <- c("HW","CW","EA","COM","SC","PER","GR","NER","CR")
abbr_legend <- c(
  "HW  = hedonic well-being",
  "CW  = cognitive well-being",
  "EA  = emotional awareness",
  "COM = compassion",
  "SC  = self-compassion",
  "PER = positive emotion regulation",
  "GR  = gratitude",
  "NER = negative emotion regulation",
  "CR  = cognitive reappraisal"
)

arm_colors <- setNames(c("#1b7837", "#762a83"), c(ARM_STR, ARM_BAR))

theme_paper <- theme_minimal(base_size = 11, base_family = FONT_FAMILY) +
  theme(text = element_text(family = FONT_FAMILY, color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, family = FONT_FAMILY),
        plot.title = element_text(size = 12, family = FONT_FAMILY),
        axis.text = element_text(size = 10, family = FONT_FAMILY, color = "black"),
        axis.title = element_text(size = 11, family = FONT_FAMILY),
        legend.text = element_text(size = 10, family = FONT_FAMILY),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA))

# ---- Figure S1: weekly means by arm ------------------------------------------
fig1 <- d %>%
  mutate(week = pmax(0, floor(days_since_start / 7))) %>%
  pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by(arm, variable, week) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            se   = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  filter(week <= 9) %>%
  mutate(variable = pretty(variable))

p1 <- ggplot(fig1, aes(week, mean, color = arm, fill = arm)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = .2, color = NA) +
  geom_line(linewidth = .7) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y") +
  scale_color_manual(values = arm_colors, name = NULL) +
  scale_fill_manual(values  = arm_colors, name = NULL) +
  scale_x_continuous(breaks = 0:9) +
  labs(x = "Week of intervention", y = "Mean (1 to 7 VAS)") +
  theme_paper

ggsave("figures/Figure_S1_trajectories.png", p1, width = 8.5, height = 7.5, dpi = 300, bg = "white")

# ---- Figure S2: arm-specific lagged slopes ------------------------------------
# Slopes and SEs from 17_phase1_slopes_by_arm.R (strength-based SE includes the covariance).
sl <- read_csv("tables/lagged_slopes_by_arm.csv", show_col_types = FALSE)
f2 <- bind_rows(
  sl %>% transmute(outcome, process, arm = ARM_BAR, beta = b_barrier,  se = se_barrier),
  sl %>% transmute(outcome, process, arm = ARM_STR, beta = b_strength, se = se_strength)) %>%
  mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         pair = paste(pretty(process), "/", pretty(outcome)))

p2 <- ggplot(f2, aes(beta, fct_rev(pair), color = arm)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .25,
                 position = position_dodge(width = .55)) +
  geom_point(size = 2.2, position = position_dodge(width = .55)) +
  scale_color_manual(values = arm_colors, name = NULL) +
  labs(x = expression(beta~"X"["t-1"]~"to"~"Y"["t+1"]~"(within-person, baseline-adjusted)"),
       y = NULL) +
  theme_paper +
  theme(axis.text.y = element_text(size = 8))

ggsave("figures/Figure_S2_lagged_effects.png", p2, width = 8.5, height = 8, dpi = 300, bg = "white")

# ---- Figure 2: temporal networks ---------------------------------------------
# Significant edges only (nonsig = "hide"), drawn when |beta| >= 0.05 (qgraph minimum).
m_cbwt <- readRDS("tables/mlvar_cbwt.rds")
m_cbt <- readRDS("tables/mlvar_cbt.rds")

mat_cbwt_t <- getNet(m_cbwt, "temporal", nonsig = "hide")
mat_cbt_t <- getNet(m_cbt, "temporal", nonsig = "hide")

rownames(mat_cbwt_t) <- colnames(mat_cbwt_t) <- abbr
rownames(mat_cbt_t) <- colnames(mat_cbt_t) <- abbr

draw_qgraph <- function(M, title_str) {
  qgraph(M,
         labels       = abbr,
         label.cex    = 1.4,
         label.scale  = FALSE,
         label.font   = 1,
         label.color  = "black",
         layout       = "circle",
         theme        = "Borkulo",
         negDashed    = FALSE,
         minimum      = 0.05,
         vsize        = 12,
         esize        = 5,
         asize        = 4,
         title        = title_str,
         title.cex    = 1.4,
         mar          = c(4, 4, 5, 4))
}

write_legend <- function() {
  par(mar = c(0.5, 0.5, 0.5, 0.5), family = FONT_FAMILY)
  plot.new()
  text(0.5, 0.97, "Node abbreviations", adj = c(0.5, 1), cex = 1.05, font = 2)
  col1 <- abbr_legend[1:5]
  col2 <- abbr_legend[6:9]
  text(0.05, 0.85, labels = paste(col1, collapse = "\n"),
       adj = c(0, 1), cex = 0.95, family = "mono")
  text(0.55, 0.85, labels = paste(col2, collapse = "\n"),
       adj = c(0, 1), cex = 0.95, family = "mono")
}

png("figures/Figure_2_mlvar_temporal.png",
    width = 3200, height = 2100, res = 320, bg = "white")
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(3, 1.0))
par(family = FONT_FAMILY)
draw_qgraph(mat_cbwt_t, "Strength-based arm")
draw_qgraph(mat_cbt_t, "Barrier-based arm")
write_legend()
dev.off()

# ---- Figure S3: change-points of the strength-based participant with most gains ----
cps <- read_csv("tables/changepoints_all.csv", show_col_types = FALSE)
ex <- cps %>% filter(is_gain, arm == "CBWT") %>%
  count(code, sort = TRUE) %>% slice(1) %>% pull(code)

d_ex <- d %>% filter(code == ex) %>% arrange(date) %>%
  mutate(t = row_number()) %>%
  select(t, all_of(c(vars, "global_wellbeing"))) %>%   # the composite is also analyzed
  pivot_longer(-t, names_to = "variable", values_to = "value") %>%
  mutate(variable = pretty(variable))

cps_ex <- cps %>% filter(code == ex, is_gain) %>%
  mutate(variable = pretty(variable))

p4 <- ggplot(d_ex, aes(t, value)) +
  geom_line(color = "grey25", linewidth = .35) +
  geom_smooth(method = "loess", se = FALSE, color = "#1b7837", span = .35, linewidth = .6) +
  geom_vline(data = cps_ex, aes(xintercept = cp_index),
             color = "#762a83", linetype = 2, linewidth = .4) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y") +
  labs(x = "Daily observation in the intervention phase", y = "Value (1 to 7)") +
  theme_paper

ggsave("figures/Figure_S3_changepoint_example.png", p4, width = 8.5, height = 7.5, dpi = 300, bg = "white")

# ---- Figure S4: participant-specific coefficients by cluster -------------------
idio   <- read_csv("tables/idionomic_betas.csv", show_col_types = FALSE)
clust  <- read_csv("tables/idionomic_clusters.csv", show_col_types = FALSE)

idio2 <- idio %>%
  inner_join(clust %>% select(code, cluster), by = "code") %>%
  mutate(arm = ifelse(arm == "CBT", ARM_BAR, ARM_STR),
         pair = paste(pretty(process), "/", pretty(outcome)),
         cluster = factor(cluster))

p5 <- ggplot(idio2, aes(pair, beta, group = code, color = cluster)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = 2) +
  geom_line(alpha = .35, linewidth = .35) +
  facet_wrap(~ arm) +
  scale_color_brewer(palette = "Dark2", name = "Cluster") +
  coord_flip() +
  labs(x = NULL, y = expression(beta~"individual lagged effect")) +
  theme_paper +
  theme(axis.text.y = element_text(size = 6.5))

ggsave("figures/Figure_S4_idionomic_clusters.png", p5, width = 8.5, height = 8, dpi = 300, bg = "white")

cat("\nFigures written to figures/.\n")
