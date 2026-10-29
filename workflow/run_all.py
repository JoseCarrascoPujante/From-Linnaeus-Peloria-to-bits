#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys, shutil
ROOT=Path(__file__).resolve().parents[1]; PY=sys.executable

def run(cmd,cwd=None):
    print('$',' '.join(map(str,cmd)),flush=True); subprocess.run([str(x) for x in cmd],check=True,cwd=cwd)

# 1 Unit test of epibit.py
epibit_root=ROOT/'01_EpiBit_Code'
run([PY,epibit_root/'tests'/'test_epibit.py'])

# 2 Raw FASTA -> count matrices (reconstructed upstream stage)
pre=ROOT/'02_FASTA_to_Counts_Reconstruction'
run([PY,pre/'code'/'reconstruct_peloria.py','--root',pre])

# 3 Deposited Supplementary Data S1, S2 and Peloria entries of S7 must match the reconstruction
run([PY,ROOT/'workflow'/'build_data_packages.py'])

# 4 Adapt reconstructed matrices to exact historical EpiBit schema
hist_in=ROOT/'05_Reproduction_Outputs'/'historical_input_matrices'; hist_in.mkdir(parents=True,exist_ok=True)
for locus in ['Lcyc','LvHIRZ']:
    run([PY,ROOT/'workflow'/'prepare_historical_matrix.py','--input',pre/'outputs'/'peloria_reconstructed'/f'{locus}_cytosine_matrix_reconstructed.csv','--output',hist_in/f'{locus}_historical_schema.csv'])

# 5 EpiBit analysis with epibit.py
epibit=epibit_root/'code'/'epibit.py'; hist_out=ROOT/'05_Reproduction_Outputs'/'historical_epibit'
for locus in ['Lcyc','LvHIRZ']:
    run([PY,epibit,'analyze','--matrix',hist_in/f'{locus}_historical_schema.csv','--outdir',hist_out/locus,'--alpha','0.5'])
run([PY,epibit,'sensitivity','--matrix',hist_in/'Lcyc_historical_schema.csv','--output',hist_out/'Lcyc_pseudocount_sensitivity.csv'])

# 6 Original epibit.py outputs versus the independent implementation
run([PY,ROOT/'workflow'/'verify_historical_epibit.py','--root',ROOT])

# 7 Comparative extension from reconstructed counts + Data S3-S7
comp=ROOT/'03_Comparative_Extension'
run([PY,comp/'code'/'calculate_epibit.py','--root',comp,'--indir',pre/'outputs'/'peloria_reconstructed','--outdir',comp/'outputs'/'epibit_metrics'])
run([PY,comp/'code'/'comparative_outputs.py','--root',comp,'--outdir',comp/'tables'])
run([PY,comp/'code'/'make_figures.py','--root',comp,'--outdir',comp/'figures'])

print('\nALL END-TO-END REPRODUCTION STEPS COMPLETED SUCCESSFULLY')
