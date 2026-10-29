#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, argparse
ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
df=pd.read_csv(a.input).rename(columns={'WT methylated':'WT_methylated','WT total':'WT_total','Peloria methylated':'Peloria_methylated','Peloria total':'Peloria_total'})
cols=['site_id','fragment','position','context','WT_methylated','WT_total','Peloria_methylated','Peloria_total']
missing=[c for c in cols if c not in df.columns]
if missing: raise SystemExit(f'Missing columns: {missing}')
Path(a.output).parent.mkdir(parents=True,exist_ok=True); df[cols].to_csv(a.output,index=False)
