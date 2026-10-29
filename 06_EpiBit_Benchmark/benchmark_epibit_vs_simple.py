#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
IN = ROOT/"inputs"
OUT = ROOT
CONTEXTS = ["CG","CHG","CHH"]
REFS = ["FWA","SUPERMAN","H06"]
NBOOT = 10000
SEED = 20261025

def prep_lcyc(df):
    a=0.5
    d=pd.DataFrame()
    d["context"]=df["context"]
    p=(df["WT methylated"]+a)/(df["WT total"]+2*a)
    q=(df["Peloria methylated"]+a)/(df["Peloria total"]+2*a)
    d["p_wt"]=p
    d["p_epi"]=q
    d["kl"]=q*np.log2(q/p)+(1-q)*np.log2((1-q)/(1-p))
    return d

def point_features(df):
    delta=df["p_epi"]-df["p_wt"]
    g=np.array([delta.mean()])
    c=np.array([delta[df["context"]==x].mean() for x in CONTEXTS])
    kl=df["kl"].to_numpy(float)
    q=np.array([kl[df["context"].to_numpy()==x].sum() for x in CONTEXTS])
    q=q/q.sum()
    p=np.array([(df["context"]==x).mean() for x in CONTEXTS])
    xy=np.array([np.log2((q[1]+q[2])/q[0]), np.log2(q[1]/p[1])])
    return {"global":g,"context":c,"epibit":xy}

def bootstrap_feature_arrays(df,nboot,seed):
    rng=np.random.default_rng(seed)
    gsum=np.zeros(nboot)
    ctx=np.zeros((nboot,3))
    kls=np.zeros((nboot,3))
    pfrac=np.zeros(3)
    nt=0
    for j,c in enumerate(CONTEXTS):
        sub=df[df["context"]==c]
        delta=(sub["p_epi"]-sub["p_wt"]).to_numpy(float)
        kl=sub["kl"].to_numpy(float)
        n=len(sub); nt+=n; pfrac[j]=n/len(df)
        idx=rng.integers(0,n,size=(nboot,n))
        ds=delta[idx]; ks=kl[idx]
        ctx[:,j]=ds.mean(axis=1)
        gsum+=ds.sum(axis=1)
        kls[:,j]=ks.sum(axis=1)
    g=(gsum/nt)[:,None]
    q=kls/kls.sum(axis=1,keepdims=True)
    xy=np.column_stack([np.log2((q[:,1]+q[:,2])/q[:,0]),np.log2(q[:,1]/pfrac[1])])
    return {"global":g,"context":ctx,"epibit":xy}

def main():
    pel=prep_lcyc(pd.read_csv(IN/"Peloria_Lcyc.csv"))
    datasets={
        "Peloria":pel,
        "FWA":pd.read_csv(IN/"FWA_merged.csv"),
        "SUPERMAN":pd.read_csv(IN/"SUPERMAN_Bondada_merged.csv"),
        "H06":pd.read_csv(IN/"H06_merged.csv"),
        "1001_Wa1":pd.read_csv(IN/"SUPERMAN_1001_Wa1.csv"),
    }
    point={k:point_features(v) for k,v in datasets.items()}
    rows=[]
    for query in ["Peloria","1001_Wa1"]:
        for method in ["global","context","epibit"]:
            q=point[query][method]
            ds={r:float(np.linalg.norm(q-point[r][method])) for r in REFS}
            rows.append({
                "query":query,"representation":method,
                "d_FWA":ds["FWA"],"d_SUPERMAN":ds["SUPERMAN"],"d_H06":ds["H06"],
                "nearest_reference":min(ds,key=ds.get),
                "H06_over_SUPERMAN_distance_ratio":ds["H06"]/ds["SUPERMAN"],
            })
    point_dist=pd.DataFrame(rows)
    point_dist.to_csv(OUT/"benchmark_point_distances.csv",index=False)

    feat=[]
    for name in datasets:
        feat.append({
            "system":name,
            "global_delta_methylation":point[name]["global"][0],
            "delta_CG":point[name]["context"][0],
            "delta_CHG":point[name]["context"][1],
            "delta_CHH":point[name]["context"][2],
            "EpiBit_X_nonCG_dominance":point[name]["epibit"][0],
            "EpiBit_Y_log2_ECHG":point[name]["epibit"][1],
        })
    pd.DataFrame(feat).to_csv(OUT/"benchmark_feature_vectors.csv",index=False)

    boot={}
    for i,(name,df) in enumerate(datasets.items()):
        boot[name]=bootstrap_feature_arrays(df,NBOOT,SEED+i)

    brows=[]; rrows=[]
    for query in ["Peloria","1001_Wa1"]:
        for method in ["global","context","epibit"]:
            q=boot[query][method]
            D=np.column_stack([np.linalg.norm(q-boot[r][method],axis=1) for r in REFS])
            winner=np.argmin(D,axis=1)
            for j,r in enumerate(REFS):
                brows.append({
                    "query":query,"representation":method,"reference":r,
                    "nearest_reference_frequency":float((winner==j).mean()),
                    "n_bootstrap":NBOOT,"seed":SEED
                })
            dsup=D[:,1]; dh06=D[:,2]; ratio=dh06/dsup
            rrows.append({
                "query":query,"representation":method,
                "P_SUPERMAN_closer_than_H06":float((dsup<dh06).mean()),
                "median_H06_over_SUPERMAN_distance_ratio":float(np.median(ratio)),
                "ratio_q05":float(np.quantile(ratio,0.05)),
                "ratio_q95":float(np.quantile(ratio,0.95)),
                "n_bootstrap":NBOOT,"seed":SEED
            })
    boot_df=pd.DataFrame(brows)
    ratio_df=pd.DataFrame(rrows)
    boot_df.to_csv(OUT/"benchmark_bootstrap_nearest_reference.csv",index=False)
    ratio_df.to_csv(OUT/"benchmark_bootstrap_distance_ratios.csv",index=False)
    point_dist.merge(
        ratio_df[["query","representation","P_SUPERMAN_closer_than_H06",
                  "median_H06_over_SUPERMAN_distance_ratio"]],
        on=["query","representation"]
    ).to_csv(OUT/"benchmark_summary.csv",index=False)

    labels={"global":"Global Δmethylation","context":"Context Δ vector","epibit":"EpiBit X/Y"}
    order=["global","context","epibit"]
    fig=plt.figure(figsize=(11,8))
    ax=fig.add_subplot(2,2,1)
    d=point_dist[point_dist["query"].eq("Peloria")].set_index("representation").loc[order]
    v=d["H06_over_SUPERMAN_distance_ratio"].to_numpy()
    ax.bar(range(3),v); ax.axhline(1,ls="--",lw=1)
    ax.set_xticks(range(3),[labels[x] for x in order],rotation=15,ha="right")
    ax.set_ylabel("H06 / SUPERMAN distance ratio")
    ax.set_title("A  Peloria: point-estimate reference separation")
    ax.text(.02,.98,">1 = SUPERMAN closer",transform=ax.transAxes,va="top",fontsize=9)

    ax=fig.add_subplot(2,2,2)
    d=ratio_df[ratio_df["query"].eq("Peloria")].set_index("representation").loc[order]
    v=100*d["P_SUPERMAN_closer_than_H06"].to_numpy()
    ax.bar(range(3),v); ax.axhline(50,ls="--",lw=1)
    ax.set_xticks(range(3),[labels[x] for x in order],rotation=15,ha="right")
    ax.set_ylabel("Bootstrap frequency (%)"); ax.set_ylim(0,100)
    ax.set_title("B  Peloria: SUPERMAN closer than H06")
    for i,x in enumerate(v): ax.text(i,x+2,f"{x:.1f}%",ha="center",fontsize=9)

    ax=fig.add_subplot(2,2,3)
    d=point_dist[point_dist["query"].eq("1001_Wa1")].set_index("representation").loc[order]
    v=d["H06_over_SUPERMAN_distance_ratio"].to_numpy()
    ax.bar(range(3),v); ax.axhline(1,ls="--",lw=1)
    ax.set_xticks(range(3),[labels[x] for x in order],rotation=15,ha="right")
    ax.set_ylabel("H06 / SUPERMAN distance ratio")
    ax.set_title("C  Independent Wa-1: point-estimate separation")
    ax.text(.02,.98,">1 = Bondada SUPERMAN closer",transform=ax.transAxes,va="top",fontsize=9)

    ax=fig.add_subplot(2,2,4)
    d=ratio_df[ratio_df["query"].eq("1001_Wa1")].set_index("representation").loc[order]
    v=100*d["P_SUPERMAN_closer_than_H06"].to_numpy()
    ax.bar(range(3),v); ax.axhline(50,ls="--",lw=1)
    ax.set_xticks(range(3),[labels[x] for x in order],rotation=15,ha="right")
    ax.set_ylabel("Bootstrap frequency (%)"); ax.set_ylim(0,100)
    ax.set_title("D  Independent Wa-1: Bondada SUPERMAN retrieval")
    for i,x in enumerate(v): ax.text(i,x+2,f"{x:.1f}%",ha="center",fontsize=9)

    fig.suptitle("EpiBit adds discriminatory information beyond conventional methylation summaries",fontsize=14)
    fig.tight_layout(rect=[0,0,1,.96])
    fig.savefig(OUT/"Figure_7_EpiBit_benchmark.png",dpi=300,bbox_inches="tight")
    plt.close(fig)

if __name__=="__main__":
    main()
