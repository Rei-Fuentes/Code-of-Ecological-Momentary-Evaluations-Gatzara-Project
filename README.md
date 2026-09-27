# Daily mechanisms of change in two interventions for chronic disease: analysis code

This repository contains every line of code behind the article *Daily Mechanisms of Change in Strength-Based and Barrier-Based Interventions for Chronic Disease: A Multilevel and Idionomic Analysis of Ecological Momentary Assessment*.

We share it so that anyone can see exactly how each number in the paper was obtained, rerun the analysis, and reuse the methods. Starting from the raw daily-assessment export, the code reproduces every number, table value and data figure in the article and its supplemental material. A verification script then recalculates 154 of those numbers and checks them against what the paper reports.

## The study in brief

In the GATZARA randomized trial, 113 adults with chronic disease took part in one of two nine-week group programs:

- **CBWT**, a contemplative, strength-based well-being training.
- **CBT**, a cognitive-behavioral, barrier-based program.

Every evening they rated nine aspects of their day on a 1–7 scale. Four of them are outcomes and six are candidate processes (global well-being is built from two of the items):

- **Outcomes:** hedonic, cognitive and global well-being, and the regulation of negative emotion.
- **Candidate processes:** emotional awareness, compassion, self-compassion, savoring, gratitude and cognitive reappraisal.

We asked how these daily processes relate to next-day well-being in each program, and we looked at the question from four angles:

| Phase | Question | Method |
|---|---|---|
| 1 | Do the day-to-day links between processes and outcomes differ between programs? | Lagged multilevel models with an arm-by-process interaction |
| 2 | How do daily states propagate through the system in each program? | Multilevel vector autoregression (temporal networks) |
| 3 | When do people change, and do processes change before outcomes? | Change-point detection (E-divisive) |
| 4 | For whom do the effects hold? | Person-specific coefficients, meta-analysis and clustering |

## What you need

- **R 4.5.1** with the packages and versions in `R_packages.txt`.
- **Python 3.12** with the packages in `requirements.txt`.

To install them:

```r
install.packages(c("dplyr", "tidyr", "readr", "lubridate", "lme4", "lmerTest", "MuMIn", "mlVAR",
                   "qgraph", "ecp", "metafor", "cluster", "simr", "ggplot2", "forcats", "stringr"))
```

```bash
pip install -r requirements.txt
```

## Data

The participant data are **not** part of this repository. The daily ratings come from people living with chronic disease, and the trial records contain personal and health information. The data are available from the corresponding author upon reasonable request.

`data/README.md` describes every file the code expects and its columns. Once you have them:

```
data/
├── raw/ema_master_completo.csv          daily EMA export (required)
└── private/
    ├── allocation.xlsx                  demographics (optional)
    ├── sample_conditions.csv            chronic-condition categories (optional)
    ├── attendance_ge4sessions.csv       session attendance (optional)
    └── code_aliases.csv                 IDs spelled differently across files (optional)
```

The `.gitignore` keeps everything in `data/raw/`, `data/private/` and `data/processed/` out of version control, as well as all generated tables, figures and logs. To add a second safeguard, run `bash tools/install_git_hook.sh` once. It installs a hook that refuses any commit containing data files.

## Running the analysis

From the repository root:

```bash
bash run_all.sh              # the whole analysis, in order
python3 verify/verify.py     # checks the reported numbers against the outputs
```

`run_all.sh` writes results to `tables/`, figures to `figures/` and one log per script to `logs/`. Every script can also be run on its own. Each one sets the random seed (20260519) at its start, so running a script alone gives the same result as running the whole pipeline.

Three steps take longer than the rest:

| Step | What it does | Time |
|---|---|---|
| `22_phase2_bootstrap.R` | 2 × 1,000 refits of the network model | about 7 h on one core, about 1 h on eight |
| `18_phase1_power_simr.R` | 5 × 1,000 simulated datasets | about 30 min |
| `33_phase3_minsize_sensitivity.R` | change-points with three window sizes | about 15 min |

Two settings control these steps:

- **`N_CORES`** spreads the bootstrap over several cores, for example `N_CORES=8 bash run_all.sh`. All resamples are drawn before any model is fitted, and the network model uses no random numbers, so the result is identical whatever the number of cores.
- **`SKIP_LONG=1`** skips the three steps and reuses their outputs if they are already in `tables/`.

## What each script does

Every script opens with a short header that explains, in plain words, what it does and why, which files it reads and writes, and where its results appear in the paper.

| Script | Purpose | Main output | In the paper |
|---|---|---|---|
| `01_preprocess.R` | Builds the analytic dataset from the raw export | `data/processed/ema_analysis_ready.csv` | Method |
| `02_sample_descriptives.py` | Sample sizes, prompts, assessment window, demographics, conditions | `sample_descriptives.csv` | Abstract, Method, Table 1, Tables S10 and S12 |
| `03_baseline_balance.R` | Did the arms start from the same place? | `baseline_balance.csv` | Results, Participant flow |
| `04_psychometrics.R` | Reliability of the single daily items | `psychometrics_*.csv` | Method, Supplementary S2 |
| `05_careless_responding.py` | Screening for careless answers | `careless_responses_*.csv` | Method, Supplementary S3, Table S9 |
| `11_phase1_lagged_mlm.R` | Phase 1: arm-by-process interactions | `phase1_lagged_mlm.csv` | Results, Phase 1, Table S1 |
| `12_phase1_sensitivity_specification.R` | Phase 1 without the baseline term, centered and standardized | `sensitivity_specification.csv` | Results, Phase 1 |
| `13_phase1_sensitivity_min30obs.R` | Phase 1 without short series | `sensitivity_min30obs.csv` | Table S7 |
| `14_phase1_sensitivity_attendance.R` | Phase 1 among participants who attended at least four sessions | `sensitivity_attendance.csv` | Table S8 |
| `15_phase1_sensitivity_careless.R` | Phase 1 without careless answers | `sensitivity_careless.csv` | Results, Phase 1, Table S9 |
| `16_phase1_sensitivity_session_group.R` | Does the session group matter? | `sensitivity_group_*.csv` | Table S14 |
| `17_phase1_slopes_by_arm.R` | Lagged slopes within each arm | `lagged_slopes_by_arm.csv` | Figure S2 |
| `18_phase1_power_simr.R` | How large an interaction could we detect? | `sensitivity_power_simr.csv` | Method, Supplementary S4 |
| `21_phase2_mlvar.R` | Phase 2: daily networks in each arm | `mlvar_*.rds`, `mlvar_*.csv` | Results, Phase 2, Figure 2 |
| `22_phase2_bootstrap.R` | Phase 2: stability of the networks | `mlvar_*_B1000*.csv` | Tables S2, S3 and S6 |
| `23_phase2_network_summaries.R` | Phase 2: density, global strength, rank stability | `network_*.csv`, `mlvar_rank_stability.csv` | Results, Phase 2, Tables S3 and S6 |
| `31_phase3_changepoints.R` | Phase 3: change-points and gains | `changepoints_*.csv` | Results, Phase 3, Table S4, Figure S3 |
| `32_phase3_gain_rate_contrast.R` | Phase 3: are gains more frequent in one arm? | `changepoint_rate_contrast.csv` | Abstract, Results, Phase 3 |
| `33_phase3_minsize_sensitivity.R` | Phase 3: sensitivity to the window size | `changepoint_minsize_sensitivity.csv` | Supplementary S7 |
| `41_phase4_idionomic.R` | Phase 4: one coefficient per person, pooled meta-analysis, clusters | `idionomic_*.csv` | Results, Phase 4, Tables S5 and S13 |
| `42_phase4_meta_by_arm.R` | Phase 4: pooled effects within each arm | `idionomic_meta_by_arm.csv` | Abstract, Results, Phase 4, Table 2 |
| `43_phase4_cluster_stability.R` | Phase 4: stability of the clusters | `cluster_bootstrap_silhouette.csv` | Supplementary S8 |
| `51_figures.R` | Figures 2 and S1–S4 | `figures/*.png` | Figures |
| `52_figure1_consort.py` | Figure 1, participant flow | `figures/Figure_1_CONSORT.png` | Figure 1 |

Outputs listed without a folder are written to `tables/`.

## Checking the results

`verify/expected_values.csv` lists 154 numbers as they appear in the manuscript and the supplement. For each one it gives the output file and the expression that recalculates it. `verify/verify.py` evaluates them all and reports PASS or FAIL. After a full run it should end with `154 passed, 0 failed`. Without the optional private files, the demographic checks are reported as SKIPPED.

## What the code does not produce

- **Formatted tables.** The tables in the article were formatted from the CSV files above and contain the same values.
- **Figure 3.** It is a conceptual diagram, not a data figure.
- **Registry counts.** The screening and allocation counts in Figure 1 come from the trial registry and are written in `52_figure1_consort.py`.
- **Randomization.** Random allocation was generated with the R package randomizr before the study started.

## Names used in the data

| Name | Meaning |
|---|---|
| `CBWT` | Strength-based arm (Contemplative Practice-Based Well-Being Training) |
| `CBT` | Barrier-based arm (cognitive-behavioral therapy) |
| `baseline`, `intervention` | Study phases; the follow-up phase is excluded |
| `session_group` | Arm and weekday of the group sessions |
| `code` | Pseudonymous participant ID |

The raw export and the trial workbooks keep their original Spanish column names. Only `01_preprocess.R` and `02_sample_descriptives.py` read them, and both translate them straight away. In the raw data the arms are labeled `EBC` and `TCC`, the Spanish acronyms of the two programs.

## Reusing this code

The code is released under the MIT license (`LICENSE`), so you may reuse and adapt it with attribution. If it helps your work, please cite the article. `session_info.txt` records the R session used for the published results.

Questions and suggestions are welcome through the issue tracker of this repository.
