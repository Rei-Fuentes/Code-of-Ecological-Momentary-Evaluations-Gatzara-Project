"""05_careless_responding.py · Did anyone answer on autopilot?

What we do here
    Daily questionnaires invite quick, careless answers. We screen every prompt with two
    standard indices (Curran, 2016; Meade & Craig, 2012): the longstring (the longest run of
    identical answers across the nine items, in item order) and the squared Mahalanobis distance
    over the same items (chi-square p value, df = 9). A prompt is flagged when longstring >= 8
    or p < .001. We exclude nothing here; 15_phase1_sensitivity_careless.R refits Phase 1
    without the prompts with longstring >= 8.

Input   data/processed/ema_analysis_ready.csv
Output  tables/careless_responses_per_prompt.csv   one row per prompt
        tables/careless_responses_summary.csv      one row per participant
Reported in  Method, Data quality; Supplementary S3; Table S9
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
SRC  = ROOT / "data" / "processed" / "ema_analysis_ready.csv"
OUT_DIR = ROOT / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# The nine EMA items
EMA_ITEMS = [
    "hedonic_wellbeing", "cognitive_wellbeing", "emotional_awareness",
    "compassion", "self_compassion", "positive_emotion_regulation",
    "gratitude", "negative_emotion_raw", "cognitive_reappraisal",
]

df = pd.read_csv(SRC)
df = df.dropna(subset=EMA_ITEMS).copy()
print(f"Prompts with all nine items: {len(df)} | participants: {df['code'].nunique()}")

# --- Longstring: longest run of identical values across the nine items ---
def longstring(row):
    vals = row[EMA_ITEMS].values
    max_run = 1; cur = 1
    for i in range(1, len(vals)):
        if vals[i] == vals[i-1]: cur += 1
        else: max_run = max(max_run, cur); cur = 1
    return max(max_run, cur)

df["longstring"] = df.apply(longstring, axis=1)

# --- Mahalanobis distance over the nine items ---
X = df[EMA_ITEMS].to_numpy()
mu = X.mean(axis=0)
cov = np.cov(X.T)
inv = np.linalg.pinv(cov)
diff = X - mu
maha2 = np.einsum("ij,jk,ik->i", diff, inv, diff)  # squared distance
df["mahalanobis_sq"] = maha2
df["p_mahalanobis"] = 1 - chi2.cdf(maha2, df=len(EMA_ITEMS))

# --- Flags per prompt ---
df["flag_longstring"] = (df["longstring"] >= 8).astype(int)  # 8 or 9 identical items in a row
df["flag_mahalanobis"] = (df["p_mahalanobis"] < .001).astype(int)
df["flag_any"] = ((df["flag_longstring"] | df["flag_mahalanobis"])).astype(int)

# --- Per-participant summary ---
summary = (df.groupby("code")
             .agg(n_prompts=("longstring", "size"),
                  longstring_mean=("longstring", "mean"),
                  longstring_max=("longstring", "max"),
                  mahalanobis_mean=("mahalanobis_sq", "mean"),
                  n_flag_longstring=("flag_longstring", "sum"),
                  n_flag_mahalanobis=("flag_mahalanobis", "sum"),
                  n_flag_any=("flag_any", "sum"))
             .reset_index()
             .sort_values("n_flag_any", ascending=False))

df.drop(columns=EMA_ITEMS).to_csv(OUT_DIR / "careless_responses_per_prompt.csv", index=False)
summary.to_csv(OUT_DIR / "careless_responses_summary.csv", index=False)

# --- Report ---
print("\n=== Careless responding ===")
print(f"Total prompts analyzed: {len(df)}")
print(f"Total participants:       {df['code'].nunique()}\n")
print("Longstring distribution (longest run of identical values across the nine items):")
print(df["longstring"].value_counts().sort_index().to_string())
print(f"\nPrompts with longstring >= 8:                  {df['flag_longstring'].sum()} "
      f"({100*df['flag_longstring'].mean():.2f}%)")
print(f"Prompts with Mahalanobis p < .001:              {df['flag_mahalanobis'].sum()} "
      f"({100*df['flag_mahalanobis'].mean():.2f}%)")
print(f"Prompts with any flag:                          {df['flag_any'].sum()} "
      f"({100*df['flag_any'].mean():.2f}%)")
print(f"\nParticipants with >=1 flagged prompt: "
      f"{(summary['n_flag_any']>0).sum()} / {len(summary)}")
print(f"Participants with >=3 flagged prompts: "
      f"{(summary['n_flag_any']>=3).sum()} / {len(summary)}")
print(f"Maximum flag rate in one participant:  "
      f"{100*(summary['n_flag_any']/summary['n_prompts']).max():.1f}%")
print("\nWritten: tables/careless_responses_per_prompt.csv, tables/careless_responses_summary.csv")
