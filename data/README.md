# Input data

The participant data are not distributed with this repository, because they come from people living
with chronic disease and include personal and health information. They are available from the
corresponding author upon reasonable request.

This page describes what the code expects, so that anyone who obtains the data, or wants to adapt
the code to their own study, knows exactly what goes where. Files placed in `data/raw/` and
`data/private/` are ignored by git.

## data/raw/ema_master_completo.csv (required)

The daily ecological momentary assessment export from the m-Path app, one row per answered prompt
and all study phases. Before this analysis it was merged with the session calendar of each group.
Its column names are the original Spanish ones; `01_preprocess.R` reads them and translates them.

| Column | Content | Becomes |
|---|---|---|
| `code` | pseudonymous participant ID | `code` |
| `tipo_intervencion` | arm: `EBC` (strength-based) or `TCC` (barrier-based) | `arm`: `CBWT` or `CBT` |
| `GRUPO` | arm and weekday of the sessions, e.g. `EBC LUNES` | `session_group`, e.g. `CBWT Monday` |
| `FECHA` | date and time the prompt was sent, `m/d/yy H:MM` | `date` |
| `fase` | `PRE`, `DURANTE` or `POST` | `phase`: `baseline`, `intervention` (follow-up dropped) |
| `dias_desde_inicio` | days from the first session of the participant's group | `days_since_start` |
| `dia_estudio` | day of study, counted from the participant's first prompt | `study_day` |
| `semana_intervencion` | intervention week, 1–9 | `intervention_week` |
| `EMA1_sliderNeutralPos` … `EMA9_sliderNeutralPos` | the nine daily items, 1–7 | named items |
| `EMA10_multipleChoice_string` | practice-time answer (not analyzed) | `practice_label` |

Phases were assigned from the days elapsed since the first group session (10, 11 or 13 February
2025 for the Monday, Tuesday and Thursday groups): baseline before it, intervention from day 0 to
day 63 (nine weekly sessions), follow-up after day 63.

## data/private/ (optional)

These files are only needed for the demographic description and for one sensitivity analysis.

**allocation.xlsx.** The trial allocation workbook, one sheet per session group (`LUNES EBC`,
`MARTES EBC`, `JUEVES EBC`, `LUNES TCC`, `MARTES TCC`, `JUEVES TCC`). `02_sample_descriptives.py`
reads four columns: `CÓDIGO` (participant ID), `EDAD` (age), `SEXO` (sex; `Mujer` = female) and
`TRASTORNO MENTAL` (comorbid mental health condition; answers starting with *s* or *y* = yes). The
script prints aggregates only.

**sample_conditions.csv.** Columns `code` and `condition_codes`: the chronic-condition categories
coded from the screening questionnaire, separated by `;`, primary condition first. Codes:
1 cardiovascular, 2 respiratory, 3 endocrine-metabolic, 4 neurological, 5 autoimmune or rheumatic,
6 gastrointestinal, 7 renal, 8 oncological, 9 chronic pain, 10 other.

**attendance_ge4sessions.csv.** One column, `code`: the participants who attended at least four of
the nine group sessions according to the trial attendance registry.

**code_aliases.csv.** Columns `workbook_code` and `ema_code`: participant IDs spelled differently in
the allocation workbook and in the EMA export. Without it one participant is not matched and the
demographics are computed on 112 of the 113 participants.
