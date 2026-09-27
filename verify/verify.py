"""verify.py — Recomputes every number of expected_values.csv from the outputs in tables/
and compares it with the value reported in the manuscript or supplement.

Run from the repository root after run_all.sh:  python3 verify/verify.py
Exit code 1 if any check fails. Checks whose inputs are missing (e.g. private files not
available) are reported as SKIPPED.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_cache = {}


def T(name):
    if name not in _cache:
        _cache[name] = pd.read_csv(ROOT / "tables" / name)
    return _cache[name]


def D(q):
    d = T("sample_descriptives.csv").set_index("quantity")["value"]
    v = d[q]
    try:
        return float(v)
    except ValueError:
        return v


def row(name, **kw):
    t = T(name)
    for k, v in kw.items():
        t = t[t[k] == v]
    assert len(t) == 1, f"{name} {kw}: {len(t)} rows"
    return t.iloc[0]


def fmt(v, dec):
    if isinstance(v, str):
        return v
    return f"{v:.{int(dec)}f}"


exp = pd.read_csv(ROOT / "verify" / "expected_values.csv", dtype=str, keep_default_na=False)
ok = fail = skip = 0
lines = []
for r in exp.itertuples():
    try:
        got = eval(r.expression, {"T": T, "D": D, "row": row, "np": np, "pd": pd})
    except (FileNotFoundError, KeyError) as e:
        skip += 1; lines.append(f"SKIP  {r.id:<6} {r.location:<32} {r.claim}  (missing: {Path(str(e.filename)).name if isinstance(e, FileNotFoundError) else e})"); continue
    got_s = fmt(got, r.decimals) if r.decimals != "" else str(got)
    passed = got_s == r.reported
    ok += passed; fail += not passed
    lines.append(f"{'PASS' if passed else 'FAIL'}  {r.id:<6} {r.location:<32} {r.claim}: reported {r.reported}, recomputed {got_s}")
print("\n".join(lines))
print(f"\n{ok} passed, {fail} failed, {skip} skipped")
sys.exit(1 if fail else 0)
