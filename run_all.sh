#!/usr/bin/env bash
# Runs the whole analysis in order, from the raw export to the figures.
#   bash run_all.sh                  everything
#   N_CORES=8 bash run_all.sh        the network bootstrap on 8 cores (same results, faster)
#   SKIP_LONG=1 bash run_all.sh      skip the three long steps and reuse their outputs in tables/
# Each script writes its log to logs/. Afterwards, check the results with:
#   python3 verify/verify.py
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p tables figures logs

run() {  # run <script>
  local s="scripts/$1"; local log="logs/${1%.*}.log"
  echo "[$(date '+%H:%M:%S')] $1"
  case "$1" in
    *.R)  Rscript "$s" > "$log" 2>&1 ;;
    *.py) python3 "$s" > "$log" 2>&1 ;;
  esac
}

run 01_preprocess.R
run 02_sample_descriptives.py
run 03_baseline_balance.R
run 04_psychometrics.R
run 05_careless_responding.py

run 11_phase1_lagged_mlm.R
run 12_phase1_sensitivity_specification.R
run 13_phase1_sensitivity_min30obs.R
run 14_phase1_sensitivity_attendance.R
run 15_phase1_sensitivity_careless.R
run 16_phase1_sensitivity_session_group.R
run 17_phase1_slopes_by_arm.R
[ "${SKIP_LONG:-0}" = 1 ] || run 18_phase1_power_simr.R            # ~30 min

run 21_phase2_mlvar.R
[ "${SKIP_LONG:-0}" = 1 ] || run 22_phase2_bootstrap.R             # several hours
run 23_phase2_network_summaries.R

run 31_phase3_changepoints.R
run 32_phase3_gain_rate_contrast.R
[ "${SKIP_LONG:-0}" = 1 ] || run 33_phase3_minsize_sensitivity.R   # ~15 min

run 41_phase4_idionomic.R
run 42_phase4_meta_by_arm.R
run 43_phase4_cluster_stability.R

run 51_figures.R
run 52_figure1_consort.py

Rscript -e 'sessionInfo()' > logs/sessionInfo.txt 2>&1
echo "Done. Outputs in tables/ and figures/, logs in logs/."
