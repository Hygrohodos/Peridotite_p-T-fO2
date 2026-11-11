#!/usr/bin/env python3
"""Compute the anchored WZ07 match level (revision v2, option B) and write data/derived/level_anchor.json."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
import pandas as pd
import make_figures as mf
from calculations.pt_fo2.level_anchor import anchor_level
from calculations.pt_fo2.params import PARAMS

df = pd.read_csv(mf.FULL, encoding="utf-8-sig", low_memory=False)
df = df.rename(columns={df.columns[0]: "SampleID"})
a = anchor_level(df)
a["params_dt_match_c"] = PARAMS.fe3.dt_match_c
(ROOT / "data" / "derived" / "level_anchor.json").write_text(json.dumps(a, indent=1))
print(json.dumps({k: v for k, v in a.items() if k != "sample_ids"}, indent=1))
