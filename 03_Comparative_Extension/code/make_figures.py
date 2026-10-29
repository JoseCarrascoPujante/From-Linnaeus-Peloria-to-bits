#!/usr/bin/env python3
from pathlib import Path
import argparse, io, zipfile
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_member(zip_path,suffix):
    with zipfile.ZipFile(zip_path) as z:
        name=[n for n in z.namelist() if n.endswith(suffix)][0]
        return pd.read_csv(io.BytesIO(z.read(name)))

def save(path):
    plt.tight_layout(); plt.savefig(path,dpi=300,bbox_inches='tight'); plt.close()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1])); ap.add_argument('--outdir',default=None); args=ap.parse_args()
    root=Path(args.root); out=Path(args.outdir or root/'figures'); out.mkdir(parents=True,exist_ok=True)
    lsite=pd.read_csv(root/'outputs'/'epibit_metrics'/'Lcyc_site_metrics.csv'); lsum=pd.read_csv(root/'outputs'/'epibit_metrics'/'Lcyc_context_summary.csv'); lvsum=pd.read_csv(root/'outputs'/'epibit_metrics'/'LvHIRZ_context_summary.csv')
    fwa=read_member(root.parent/'04_Supplementary_Data'/'Data_S3_FWA.zip','FWA_summary.csv'); sup=read_member(root.parent/'04_Supplementary_Data'/'Data_S4_SUPERMAN.zip','SUPERMAN_summary.csv'); h06=read_member(root.parent/'04_Supplementary_Data'/'Data_S6_H06.zip','H06_summary.csv'); h06w=read_member(root.parent/'04_Supplementary_Data'/'Data_S6_H06.zip','H06_window_sensitivity.csv'); arch=read_member(root.parent/'04_Supplementary_Data'/'Data_S7_architecture_space.zip','architecture_source_and_coordinates.csv'); perm=read_member(root.parent/'04_Supplementary_Data'/'Data_S7_architecture_space.zip','context_label_permutation_summary.csv')
    colors={'CG':'#4472C4','CHG':'#70AD47','CHH':'#ED7D31'}

    # Figure 2: data-only reconstruction of quantitative panels
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    x=np.arange(1,len(lsite)+1)
    for ctx in ['CG','CHG','CHH']:
        s=lsite[lsite.context==ctx]; axs[0].scatter(s.index+1,s['WT p'].astype(float),s=18,label=f'WT {ctx}',color=colors[ctx],alpha=.65); axs[0].scatter(s.index+1,s['Peloria q'].astype(float),s=18,marker='x',label=f'Peloria {ctx}',color=colors[ctx])
    axs[0].set(xlabel='Informative cytosine',ylabel='Methylation probability',ylim=(-.02,1.02),title='Cytosine-level methylation probabilities'); axs[0].legend(fontsize=7,ncol=2,frameon=False)
    s=lsum[lsum.context.isin(['CG','CHG','CHH'])]; xx=np.arange(3); w=.36; axs[1].bar(xx-w/2,s.WT_mean,w,label='WT'); axs[1].bar(xx+w/2,s.Peloria_mean,w,label='Peloria'); axs[1].set_xticks(xx,s.context); axs[1].set(ylabel='Mean corrected methylation',title='Context means'); axs[1].legend(frameon=False)
    fig.suptitle('Figure 2 source reconstruction: Lcyc cytosine-resolved methylation'); save(out/'Figure_2_reproduced.png')

    # Figure 3
    fig,axs=plt.subplots(1,3,figsize=(11,4)); s=lsum[lsum.context.isin(['CG','CHG','CHH'])]; xx=np.arange(3); w=.36
    axs[0].bar(xx-w/2,s.H_WT,w,label='WT'); axs[0].bar(xx+w/2,s.H_Peloria,w,label='Peloria'); axs[0].set_xticks(xx,s.context); axs[0].set(ylabel='Summed entropy (bits)',title='Shannon entropy'); axs[0].legend(frameon=False)
    total=lsum[lsum.context=='Total'].iloc[0]; axs[1].bar(['KL','JSD','Dsym'],[total.KL_Peloria_WT,total.JSD,total.D_sym]); axs[1].set(ylabel='Bits',title='Regional divergence')
    axs[2].bar(s.context,s.KL_Peloria_WT,color=[colors[c] for c in s.context]); axs[2].set(ylabel='Directional KL (bits)',title='Context contributions')
    fig.suptitle('Figure 3 reproduced from reconstructed Lcyc matrix'); save(out/'Figure_3_reproduced.png')

    # Figure 4
    ranked=lsite.sort_values('D_sym',ascending=False).reset_index(drop=True); cum=ranked.D_sym.astype(float).cumsum()/ranked.D_sym.astype(float).sum()*100
    fig,axs=plt.subplots(1,2,figsize=(9,4)); axs[0].plot(np.arange(1,len(ranked)+1),cum); axs[0].set(xlabel='Top-ranked cytosines',ylabel='Cumulative Dsym (%)',title='Cumulative divergence');
    for t in [50,80,90,95]: axs[0].axhline(t,ls='--',lw=.7,color='grey')
    comp=pd.DataFrame({'Locus':['Lcyc','LvHIRZ'],'KL':[total.KL_Peloria_WT,lvsum[lvsum.context=='Total'].iloc[0].KL_Peloria_WT],'JSD':[total.JSD,lvsum[lvsum.context=='Total'].iloc[0].JSD],'Dsym':[total.D_sym,lvsum[lvsum.context=='Total'].iloc[0].D_sym]}); comp.set_index('Locus').plot(kind='bar',ax=axs[1]); axs[1].set_yscale('log'); axs[1].set(ylabel='Bits (log scale)',title='Locus specificity'); axs[1].tick_params(axis='x',rotation=0)
    fig.suptitle('Figure 4 reproduced: dimensionality and locus specificity'); save(out/'Figure_4_reproduced.png')

    # Figure 5: pseudocount sensitivity recomputed from counts
    alphas=[.1,.25,.5,1,2]; rows=[]
    for a in alphas:
        p=(lsite['WT methylated'].astype(float)+a)/(lsite['WT total'].astype(float)+2*a); q=(lsite['Peloria methylated'].astype(float)+a)/(lsite['Peloria total'].astype(float)+2*a)
        kl=q*np.log2(q/p)+(1-q)*np.log2((1-q)/(1-p)); rev=p*np.log2(p/q)+(1-p)*np.log2((1-p)/(1-q)); ds=(kl+rev)/2; m=(p+q)/2; js=.5*(p*np.log2(p/m)+(1-p)*np.log2((1-p)/(1-m))+q*np.log2(q/m)+(1-q)*np.log2((1-q)/(1-m)))
        chg=lsite.context.eq('CHG'); noncg=~lsite.context.eq('CG'); rank=np.sort(ds.to_numpy())[::-1]; cumulative=np.cumsum(rank)/rank.sum(); nth=[int(np.argmax(cumulative>=t))+1 for t in [.5,.8,.9,.95]]
        rows.append([a,kl.sum(),ds.sum(),js.sum(),kl[chg].sum()/kl.sum()*100,kl[noncg].sum()/kl.sum()*100,*nth])
    sens=pd.DataFrame(rows,columns=['alpha','KL','Dsym','JSD','CHG_pct','nonCG_pct','N50','N80','N90','N95']); sens.to_csv(root/'tables'/'Main_Table_2_pseudocount_sensitivity.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(9,4)); axs[0].plot(sens.alpha,sens.KL,label='KL'); axs[0].plot(sens.alpha,sens.Dsym,label='Dsym'); axs[0].plot(sens.alpha,sens.JSD,label='JSD'); axs[0].set(xlabel='alpha',ylabel='Bits',title='Divergence sensitivity'); axs[0].legend(frameon=False); axs[1].plot(sens.alpha,sens.CHG_pct,label='CHG % KL'); axs[1].plot(sens.alpha,sens.nonCG_pct,label='non-CG % KL'); axs[1].set(xlabel='alpha',ylabel='Percent',title='Context contribution sensitivity'); axs[1].legend(frameon=False); fig.suptitle('Figure 5 reproduced: pseudocount robustness'); save(out/'Figure_5_reproduced.png')

    # Figure 6
    fig,axs=plt.subplots(2,2,figsize=(9,7));
    for _,r in arch.iterrows(): axs[0,0].scatter(r.X_nonCG_dominance,r.Y_log2_ECHG,s=70); axs[0,0].text(r.X_nonCG_dominance+.08,r.Y_log2_ECHG+.05,r.system,fontsize=8)
    axs[0,0].axhline(0,ls='--',color='grey',lw=.7); axs[0,0].axvline(0,ls='--',color='grey',lw=.7); axs[0,0].set(xlabel='non-CG dominance X',ylabel='CHG enrichment Y',title='Architecture space')
    arch.set_index('system')[['q_CG','q_CHG','q_CHH']].plot(kind='bar',stacked=True,ax=axs[0,1],color=[colors['CG'],colors['CHG'],colors['CHH']]); axs[0,1].set(title='Context allocation',ylabel='q_c'); axs[0,1].tick_params(axis='x',rotation=25)
    arch.set_index('system')[['E_CG','E_CHG','E_CHH']].plot(kind='bar',ax=axs[1,0],color=[colors['CG'],colors['CHG'],colors['CHH']]); axs[1,0].axhline(1,ls='--',color='grey',lw=.7); axs[1,0].set(title='Context enrichment',ylabel='E_c'); axs[1,0].tick_params(axis='x',rotation=25)
    for _,r in arch.iterrows(): axs[1,1].plot([50,80,90,95],[r.R50,r.R80,r.R90,r.R95],marker='o',label=r.system)
    axs[1,1].set(xlabel='Percent Dsym explained',ylabel='Fraction of sites required',title='Normalized dimensionality'); axs[1,1].legend(frameon=False,fontsize=8)
    fig.suptitle('Figure 6 reproduced: epiallelic architecture space'); save(out/'Figure_6_reproduced.png')

    # Supplementary data-driven figures S2-S13, source-oriented compact versions
    # S2 is Figure 2 probability panel
    # copy by regeneration to distinct names
    for src,dst in [('Figure_2_reproduced.png','Supplementary_Figure_S2_reproduced.png'),('Figure_4_reproduced.png','Supplementary_Figure_S4_reproduced.png'),('Figure_6_reproduced.png','Supplementary_Figure_S9_reproduced.png')]:
        (out/dst).write_bytes((out/src).read_bytes())
    # S3 context fractions
    fig,ax=plt.subplots(figsize=(6,4)); ss=lsum[lsum.context.isin(['CG','CHG','CHH'])]; ax.bar(ss.context,ss.KL_Peloria_WT,color=[colors[c] for c in ss.context]); ax.set(ylabel='Directional KL bits',title='Supplementary Figure S3 reproduced: Lcyc context contributions'); save(out/'Supplementary_Figure_S3_reproduced.png')
    # S5 locus specificity
    (out/'Supplementary_Figure_S5_reproduced.png').write_bytes((out/'Figure_4_reproduced.png').read_bytes())
    # S6 FWA summary
    fig,axs=plt.subplots(1,2,figsize=(8,4)); fwa.set_index('analysis')[['q_CG','q_CHG','q_CHH']].plot(kind='bar',stacked=True,ax=axs[0],color=[colors['CG'],colors['CHG'],colors['CHH']]); fwa.set_index('analysis')[['E_CG','E_CHG','E_CHH']].plot(kind='bar',ax=axs[1],color=[colors['CG'],colors['CHG'],colors['CHH']]); axs[1].axhline(1,ls='--',color='grey'); fig.suptitle('Supplementary Figure S6 reproduced: FWA'); save(out/'Supplementary_Figure_S6_reproduced.png')
    fig,axs=plt.subplots(1,2,figsize=(8,4)); sup.set_index('analysis')[['q_CG','q_CHG','q_CHH']].plot(kind='bar',stacked=True,ax=axs[0],color=[colors['CG'],colors['CHG'],colors['CHH']]); sup.set_index('analysis')[['E_CG','E_CHG','E_CHH']].plot(kind='bar',ax=axs[1],color=[colors['CG'],colors['CHG'],colors['CHH']]); axs[1].axhline(1,ls='--',color='grey'); fig.suptitle('Supplementary Figure S7/S10 reproduced: SUPERMAN'); save(out/'Supplementary_Figure_S7_reproduced.png'); (out/'Supplementary_Figure_S10_reproduced.png').write_bytes((out/'Supplementary_Figure_S7_reproduced.png').read_bytes())
    fig,axs=plt.subplots(1,2,figsize=(8,4)); h06.set_index('analysis')[['q_CG','q_CHG','q_CHH']].plot(kind='bar',stacked=True,ax=axs[0],color=[colors['CG'],colors['CHG'],colors['CHH']]); h06.set_index('analysis')[['E_CG','E_CHG','E_CHH']].plot(kind='bar',ax=axs[1],color=[colors['CG'],colors['CHG'],colors['CHH']]); axs[1].axhline(1,ls='--',color='grey'); fig.suptitle('Supplementary Figure S8 reproduced: H06'); save(out/'Supplementary_Figure_S8_reproduced.png')
    # S11 from stored exact permutation summary: observed vs null 0.1-99.9% bands
    fig,ax=plt.subplots(figsize=(9,5)); pp=perm.copy(); x=np.arange(len(pp)); ax.vlines(x,pp.null_p0_1,pp.null_p99_9,color='lightgrey',lw=7,label='0.1–99.9% null interval'); ax.scatter(x,pp.observed_divergence_fraction,c='black',s=25,zorder=3,label='Observed q_c'); ax.set_xticks(x,[f'{a}\n{b}' for a,b in zip(pp.system,pp.context)],rotation=45,ha='right'); ax.set_ylabel('Directional-KL divergence fraction'); ax.set_title('Supplementary Figure S11 reproduced: context-label permutation tests'); ax.legend(frameon=False); save(out/'Supplementary_Figure_S11_reproduced.png')
    fig,axs=plt.subplots(1,2,figsize=(9,4)); h06w.set_index('window')[['q_CG','q_CHG','q_CHH']].plot(kind='bar',stacked=True,ax=axs[0],color=[colors['CG'],colors['CHG'],colors['CHH']]); h06w.set_index('window')[['E_CG','E_CHG','E_CHH']].plot(kind='bar',ax=axs[1],color=[colors['CG'],colors['CHG'],colors['CHH']]); axs[1].axhline(1,ls='--',color='grey'); fig.suptitle('Supplementary Figure S12 reproduced: H06 window sensitivity'); save(out/'Supplementary_Figure_S12_reproduced.png')
    fig,ax=plt.subplots(figsize=(6,5));
    for _,r in arch.iterrows(): ax.scatter(r.ILR1_CG_vs_nonCG,r.ILR2_CHG_vs_CHH,s=70); ax.text(r.ILR1_CG_vs_nonCG+.05,r.ILR2_CHG_vs_CHH+.03,r.system,fontsize=8)
    ax.axhline(0,ls='--',color='grey',lw=.7); ax.axvline(0,ls='--',color='grey',lw=.7); ax.set(xlabel='ILR1 CG vs non-CG',ylabel='ILR2 CHG vs CHH',title='Supplementary Figure S13 reproduced: compositional sensitivity'); save(out/'Supplementary_Figure_S13_reproduced.png')
    print(f'Wrote reproducible data-driven figures to {out}')
if __name__=='__main__':main()
