"""02_sample_descriptives.py · Who took part, and how much data did they give?

What we do here
    We compute every descriptive number of the sample that the paper reports: participants per
    arm, prompts, prompts per participant, the assessment window, how many lags join
    consecutive days, the effective sample of each analysis, and, when the private files are
    available, demographics and chronic-condition categories.
    Only aggregates are written or printed. The private files contain personal and health data
    and never leave data/private/.

Input   data/processed/ema_analysis_ready.csv
        data/private/allocation.xlsx         (optional, not distributed) age, sex, comorbidity
        data/private/sample_conditions.csv   (optional, not distributed) condition categories
        data/private/code_aliases.csv        (optional, not distributed) IDs spelled differently
Output  tables/sample_descriptives.csv       one row per reported quantity
Reported in  Abstract; Method; Results, Participant flow; Table 1; Tables S10 and S12

Run from the repository root:  python3 scripts/02_sample_descriptives.py
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tables"
OUT.mkdir(exist_ok=True)
rows = []


def add(quantity, value, where):
    rows.append({"quantity": quantity, "value": value, "reported_in": where})


d = pd.read_csv(ROOT / "data" / "processed" / "ema_analysis_ready.csv", parse_dates=["date"])
arm_label = {"CBWT": "strength_based", "CBT": "barrier_based"}

# --- Analytic sample and prompts -----------------------------------------------------------
add("participants", d.code.nunique(), "Method; Results, Participant flow")
for arm, lab in arm_label.items():
    add(f"participants_{lab}", d.loc[d.arm == arm, "code"].nunique(), "Method; Results")
add("daily_prompts_total", len(d), "Method, Data quality; Results; Table S9")
add("baseline_prompts", int((d.phase == "baseline").sum()), "Supplementary")

dur = d[d.phase == "intervention"]
n_dur = dur.groupby(["arm", "code"]).size()
add("intervention_obs_median_all", float(n_dur.median()), "Abstract; Supplementary S2, S4")
for arm, lab in arm_label.items():
    s = n_dur.loc[arm]
    add(f"intervention_obs_{lab}", int(s.sum()), "Table 1")
    add(f"intervention_obs_median_{lab}", float(s.median()), "Results, Participant flow; Table 1; Figure 1")
    add(f"intervention_obs_range_{lab}", f"{s.min()}-{s.max()}", "Table 1")

# --- Analysis-specific samples (Table S12) ------------------------------------------------
for thr, where in ((20, "Phase 2 (mlVAR); Table S12"), (28, "Phase 3 (change-points); Table S12"),
                   (30, "Phase 1 sensitivity; Table S7")):
    keep = n_dur[n_dur >= thr]
    add(f"participants_ge{thr}obs", len(keep), where)
    for arm, lab in arm_label.items():
        add(f"participants_ge{thr}obs_{lab}", int(len(keep.loc[arm])) if arm in keep.index.get_level_values(0) else 0, where)
        if thr == 20:
            add(f"intervention_obs_ge20_{lab}", int(keep.loc[arm].sum()), "Supplementary S6 (observations entering mlVAR)")

# --- Assessment window ---------------------------------------------------------------------
span = d.groupby("code").date.agg(["min", "max"])
elapsed = (span["max"] - span["min"]).dt.days
add("window_days_median", float(elapsed.median()), "Method, Measures (days elapsed, first to last prompt)")
add("window_days_range", f"{elapsed.min()}-{elapsed.max()}", "Method, Measures (days elapsed)")
for ph, lab in (("baseline", "baseline"), ("intervention", "intervention")):
    s = d[d.phase == ph].groupby("code").date.agg(["min", "max"])
    add(f"{lab}_span_days_median", float(((s["max"] - s["min"]).dt.days + 1).median()),
        "Method (calendar days, inclusive)")

# --- Lag coverage --------------------------------------------------------------------------
consec = d.consecutive_day.astype(str).str.lower().eq("true")
add("lag_coverage_pct", round(100 * consec.mean(), 1), "Method, Data processing")

# --- Demographics (private allocation workbook) --------------------------------------------
alloc = ROOT / "data" / "private" / "allocation.xlsx"
if alloc.exists():
    norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())
    # Optional: IDs written differently in the workbook and in the EMA export (columns
    # workbook_code, ema_code), kept with the private files.
    alias_f = ROOT / "data" / "private" / "code_aliases.csv"
    alias = ({norm(a): norm(b) for a, b in pd.read_csv(alias_f).itertuples(index=False)}
             if alias_f.exists() else {})
    sheets = ["LUNES EBC", "MARTES EBC", "JUEVES EBC", "LUNES TCC", "MARTES TCC", "JUEVES TCC"]
    a = pd.concat([pd.read_excel(alloc, sheet_name=s) for s in sheets])
    a["_c"] = a["CÓDIGO"].map(norm).replace(alias)
    codes = {norm(c) for c in d.code.unique()}
    s = a[a["_c"].isin(codes)].drop_duplicates("_c")
    age = pd.to_numeric(s["EDAD"], errors="coerce")
    sex = s["SEXO"].astype(str).str.strip().str.lower()
    mh = s["TRASTORNO MENTAL"].astype(str).str.strip().str.lower().str.startswith(("s", "y"))
    add("demographics_matched", len(s), "check: must equal participants")
    add("age_mean", round(age.mean(), 1), "Method, Participants")
    add("age_sd", round(age.std(), 2), "Method, Participants")
    add("age_range", f"{age.min():.0f}-{age.max():.0f}", "Method, Participants")
    add("female_pct", round(100 * (sex == "mujer").mean(), 1), "Method, Participants; Limitations")
    add("mental_health_comorbidity_pct", round(100 * mh.mean(), 1), "Method, Participants")
else:
    print("data/private/allocation.xlsx not found: demographics skipped")

# --- Chronic condition categories (Table S10) ----------------------------------------------
cond = ROOT / "data" / "private" / "sample_conditions.csv"
if cond.exists():
    cats = {1: "Cardiovascular", 2: "Respiratory", 3: "Endocrine-metabolic", 4: "Neurological",
            5: "Autoimmune / rheumatic", 6: "Gastrointestinal", 7: "Renal", 8: "Oncological",
            9: "Chronic pain", 10: "Other"}
    c = pd.read_csv(cond, dtype={"condition_codes": str})
    c = c.merge(d[["code", "arm"]].drop_duplicates(), on="code")
    # Primary condition = first code listed (participants with comorbid conditions keep all codes)
    c["primary"] = c.condition_codes.str.split(";").str[0].astype(int).map(cats)
    n = len(c)
    for cat, k in c.primary.value_counts().items():
        add(f"condition_{cat}_n", int(k), "Table S10")
        add(f"condition_{cat}_pct", round(100 * k / n, 1), "Method, Participants; Table S10")
        for arm, lab in arm_label.items():
            add(f"condition_{cat}_{lab}_n", int(((c.primary == cat) & (c.arm == arm)).sum()), "Table S10")
else:
    print("data/private/sample_conditions.csv not found: Table S10 skipped")

out = pd.DataFrame(rows)
out.to_csv(OUT / "sample_descriptives.csv", index=False)
print(out.to_string(index=False))
