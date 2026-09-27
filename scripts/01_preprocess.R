# =============================================================================
# 01_preprocess.R · From the raw EMA export to the analytic dataset
# =============================================================================
# What we do here
#   We take the daily ecological momentary assessment (EMA) export and build the one table
#   that every later script reads. The raw export keeps the original Spanish column names and
#   labels; this is the only script that sees them, and it hands everything on in English.
#
#   1. Give the nine daily items readable names.
#   2. Reverse-score the negative-emotion item, so that higher always means better.
#   3. Build global well-being as the mean of hedonic and cognitive well-being.
#   4. Keep the baseline and intervention phases (follow-up had about two prompts per person).
#   5. Keep one prompt per participant per day and flag rows that follow the previous
#      calendar day, which the lagged models need.
#   6. Keep participants with at least one intervention-phase prompt (114 -> 113).
#
# Input   data/raw/ema_master_completo.csv      (not distributed; see data/README.md)
# Output  data/processed/ema_analysis_ready.csv
#         data/processed/data_dictionary.md
#         data/processed/preprocess_log.txt
# Reported in  Method (Measures; Data quality and diagnostic procedures; Data processing)
#
# Run from the repository root:  Rscript scripts/01_preprocess.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
})

raw_path <- "data/raw/ema_master_completo.csv"
out_dir  <- "data/processed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_lines <- c()
log <- function(msg) {
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log(paste("== Preprocessing ==", Sys.time()))

# -----------------------------------------------------------------------------
# 1. Load and translate the columns we use
# -----------------------------------------------------------------------------
raw <- read_csv(raw_path, show_col_types = FALSE)
log(sprintf("Loaded: %d rows, %d columns, %d participants",
            nrow(raw), ncol(raw), n_distinct(raw$code)))

weekday <- c(LUNES = "Monday", MARTES = "Tuesday", JUEVES = "Thursday")

df <- raw %>%
  transmute(
    code,
    # Arms: EBC and TCC are the Spanish acronyms of the two programs.
    arm = recode(tipo_intervencion, EBC = "CBWT", TCC = "CBT"),
    session_group = paste(recode(sub(" .*", "", GRUPO), EBC = "CBWT", TCC = "CBT"),
                          weekday[sub(".* ", "", GRUPO)]),
    sent = FECHA,
    phase = recode(fase, PRE = "baseline", DURANTE = "intervention", POST = "follow_up"),
    study_day = dia_estudio,
    days_since_start = dias_desde_inicio,
    intervention_week = semana_intervencion,
    hedonic_wellbeing           = EMA1_sliderNeutralPos,   # happiness today
    cognitive_wellbeing         = EMA2_sliderNeutralPos,   # day satisfaction
    emotional_awareness         = EMA3_sliderNeutralPos,   # mindful awareness
    compassion                  = EMA4_sliderNeutralPos,   # compassion toward others
    self_compassion             = EMA5_sliderNeutralPos,
    positive_emotion_regulation = EMA6_sliderNeutralPos,   # savoring of positive events
    gratitude                   = EMA7_sliderNeutralPos,
    negative_emotion_raw        = EMA8_sliderNeutralPos,   # difficulty recovering from negative emotion
    cognitive_reappraisal       = EMA9_sliderNeutralPos,
    practice_label              = EMA10_multipleChoice_string
  )

# -----------------------------------------------------------------------------
# 2-3. Reverse the negative-emotion item and build the composite
# -----------------------------------------------------------------------------
df <- df %>%
  mutate(negative_emotion_regulation = 8 - negative_emotion_raw,
         global_wellbeing = (hedonic_wellbeing + cognitive_wellbeing) / 2)
log(sprintf("Negative-emotion item reversed: original mean %.2f, reversed %.2f",
            mean(df$negative_emotion_raw, na.rm = TRUE),
            mean(df$negative_emotion_regulation, na.rm = TRUE)))
log(sprintf("Global well-being: mean %.2f, sd %.2f",
            mean(df$global_wellbeing, na.rm = TRUE), sd(df$global_wellbeing, na.rm = TRUE)))

# Practice time (tenth item, not analyzed) recoded to minutes. The answer labels are Spanish.
df <- df %>%
  mutate(practice_time = case_when(
    grepl("No he practicado", practice_label) ~ 0L,   # "I did not practice"
    grepl("0-10 min",         practice_label) ~ 5L,
    grepl("10 - 20 min",      practice_label) ~ 15L,
    grepl("20 - 30 min",      practice_label) ~ 25L,
    grepl("s de 30 min",      practice_label) ~ 35L,  # "Más de 30 min", matched without the accent
    TRUE                                      ~ NA_integer_
  ))

# -----------------------------------------------------------------------------
# 4. Baseline and intervention phases only
# -----------------------------------------------------------------------------
n_before <- nrow(df)
df <- df %>% filter(phase %in% c("baseline", "intervention"))
log(sprintf("Phase filter: %d -> %d rows (baseline = %d, intervention = %d)",
            n_before, nrow(df), sum(df$phase == "baseline"), sum(df$phase == "intervention")))

# -----------------------------------------------------------------------------
# 5. One prompt per day and consecutive-day flag
# -----------------------------------------------------------------------------
df <- df %>% mutate(date = as.Date(mdy_hm(sent)))

dup <- df %>% count(code, date) %>% filter(n > 1)
log(sprintf("Participant-days with more than one prompt: %d (of %d); we keep the first.",
            nrow(dup), n_distinct(paste(df$code, df$date))))

df <- df %>%
  arrange(code, date) %>%
  group_by(code, date) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(code, date) %>%
  group_by(code) %>%
  mutate(
    days_since_previous = as.integer(date - lag(date)),
    consecutive_day     = !is.na(days_since_previous) & days_since_previous == 1L
  ) %>%
  ungroup()
log(sprintf("After deduplication: %d rows, %d participants", nrow(df), n_distinct(df$code)))

# -----------------------------------------------------------------------------
# 6. Minimum inclusion criterion: at least one intervention-phase prompt
# -----------------------------------------------------------------------------
with_intervention <- df %>% filter(phase == "intervention") %>% distinct(code) %>% pull(code)
n_drop <- n_distinct(df$code) - length(with_intervention)
df <- df %>% filter(code %in% with_intervention)
log(sprintf("Participants without intervention-phase data dropped: %d -> %d rows, %d participants",
            n_drop, nrow(df), n_distinct(df$code)))
log(sprintf("Rows that follow the previous calendar day: %.1f%%", 100 * mean(df$consecutive_day)))

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
df_out <- df %>%
  select(code, arm, session_group, date, phase, study_day, days_since_start,
         intervention_week, consecutive_day, days_since_previous,
         hedonic_wellbeing, cognitive_wellbeing, global_wellbeing,
         emotional_awareness, compassion, self_compassion,
         positive_emotion_regulation, gratitude,
         negative_emotion_raw, negative_emotion_regulation, cognitive_reappraisal,
         practice_time, practice_label)

by_arm <- df_out %>%
  count(arm, phase, code) %>%
  group_by(arm, phase) %>%
  summarise(prompts = sum(n), participants = n(), median_prompts = median(n), .groups = "drop")
log("\nBy arm and phase:")
log(paste(capture.output(print(by_arm)), collapse = "\n"))

write_csv(df_out, file.path(out_dir, "ema_analysis_ready.csv"))
log(sprintf("\nWritten: %s", file.path(out_dir, "ema_analysis_ready.csv")))
writeLines(log_lines, file.path(out_dir, "preprocess_log.txt"))

dict <- c(
  "# Data dictionary: ema_analysis_ready.csv",
  "",
  sprintf("Rows: %d  |  Participants: %d", nrow(df_out), n_distinct(df_out$code)),
  "",
  "| Column | Type | Description |",
  "|---|---|---|",
  "| code | text | Pseudonymous participant ID |",
  "| arm | text | CBWT = strength-based (Contemplative Practice-Based Well-Being Training); CBT = barrier-based (cognitive-behavioral therapy) |",
  "| session_group | text | Arm and weekday of the group sessions (Monday, Tuesday or Thursday) |",
  "| date | date | Date of the prompt |",
  "| phase | text | baseline or intervention |",
  "| study_day | integer | Day of study, counted from the participant's first prompt (1 = first) |",
  "| days_since_start | integer | Days since the first session of the participant's group (negative at baseline) |",
  "| intervention_week | integer | Intervention week, 1-9 (empty at baseline) |",
  "| consecutive_day | logical | TRUE if the participant's previous row is the previous calendar day |",
  "| days_since_previous | integer | Days since the participant's previous prompt |",
  "| hedonic_wellbeing | 1-7 | Happiness today |",
  "| cognitive_wellbeing | 1-7 | Day satisfaction |",
  "| global_wellbeing | 1-7 | Mean of hedonic and cognitive well-being |",
  "| emotional_awareness | 1-7 | Mindful awareness |",
  "| compassion | 1-7 | Compassion toward others |",
  "| self_compassion | 1-7 | Self-compassion |",
  "| positive_emotion_regulation | 1-7 | Savoring of positive events |",
  "| gratitude | 1-7 | Gratitude |",
  "| negative_emotion_raw | 1-7 | Difficulty recovering from negative emotion, as answered |",
  "| negative_emotion_regulation | 1-7 | 8 minus the previous column (higher = better regulation) |",
  "| cognitive_reappraisal | 1-7 | Cognitive reappraisal |",
  "| practice_time | integer | Minutes of practice that day (0, 5, 15, 25 or 35); not analyzed |",
  "| practice_label | text | Original answer label of the practice item (Spanish) |"
)
writeLines(dict, file.path(out_dir, "data_dictionary.md"))
log("\n== Preprocessing finished ==")
