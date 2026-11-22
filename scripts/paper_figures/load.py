"""Repo access for the paper-figure scripts. REPO defaults to this checkout; set PAPER_FIG_REPO to point elsewhere."""
import sys, os
from pathlib import Path
sys.dont_write_bytecode=True
REPO=os.environ.get('PAPER_FIG_REPO', str(Path(__file__).resolve().parents[2]))
sys.path.insert(0,REPO); sys.path.insert(0,REPO+'/scripts')
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import make_figures as mf          # helper functions only; nothing is written to the repo
from calculations.pt_fo2.params import PARAMS
from calculations.pt_fo2 import calibration as calib
def load():
    df=pd.read_csv(mf.FULL,encoding='utf-8-sig',low_memory=False)
    return df.rename(columns={df.columns[0]:'SampleID'})
