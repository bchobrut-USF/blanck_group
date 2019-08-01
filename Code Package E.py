import pandas as pd
from scipy.stats.stats import pearsonr
from lifelines import CoxPHFitter
from multiprocessing import Pool, cpu_count
import pandas as pd
from scipy import stats
import numpy as np
import sys
import os
from lifelines.statistics import logrank_test
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu

def coxcalc(df, x, survivaltime, status):
    df5 = df[[status, survivaltime, x]]
    df5[x] = pd.to_numeric(df5[x])
    df5 = df5.dropna()
    cph = CoxPHFitter()
    cph.fit(df5, duration_col=survivaltime, event_col=status, show_progress=False)
    return cph.summary


def get_data():
    rna = pd.read_csv("")
    rna = rna.drop_duplicates("Hugo_Symbol")
    rna = rna.transpose()
    rna.columns = rna.iloc[0]
    rna = rna.reindex(rna.index.drop("Hugo_Symbol"))
    rna.index = rna.index.str.slice(0,12)
    rna = rna[~rna.index.duplicated(keep='first')]
    
    comp = pd.read_csv("")
    comp = comp[comp["Hugo_Symbol"] == "IDH1"]
    comp = comp.set_index("Filename")
    comp = comp[["CSU", "CSN", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]]
    
    df = comp.join(rna)
    
    return df


def mp_calc(gene):
    print gene
    global df

    try:
        df1 = df[df["CSN"] > 0]
        df2 = df[df["CSN"] <= 0]
        df3 = df[df["CSU"] >= df["CSU"].median()]
        df4 = df[df["CSU"] < df["CSU"].median()]
        mwu_ncpr = mannwhitneyu(df1[gene], df2[gene])
        mwu_uversky = mannwhitneyu(df3[gene], df4[gene])
        
        ncpr_comp_rna = df1[gene].mean()
        ncpr_noncomp_rna = df2[gene].mean()
        
        uv_comp_rna = df3[gene].mean()
        uv_noncomp_rna = df4[gene].mean()
        
        return gene, mwu_ncpr[1], mwu_uversky[1], ncpr_comp_rna, ncpr_noncomp_rna, uv_comp_rna, uv_noncomp_rna
    except:
        print 'error'
        result = pd.DataFrame()
        return None, None, None, None, None, None, None
gene_list = []
ncpr_ps = []
uversky_ps = []
ncpr_comp_rna_list = []
ncpr_noncomp_rna_list = []
uv_comp_rna_list = []
uv_noncomp_rna_list = []

def rna_comp_corr():
    global df
    genes = list(df)
    genes.remove("CSU")
    genes.remove("CSN")
    genes.remove("OS_MONTHS")
    genes.remove("OS_STATUS")
    genes.remove("DFS_MONTHS")
    genes.remove("DFS_STATUS")
    
    results = pd.DataFrame()
    
    p = Pool(cpu_count())
    source = p.map(mp_calc, genes)
    for gene, logrank_ncpr, logrank_uversky, ncpr_comp_rna, ncpr_noncomp_rna, uv_comp_rna, uv_noncomp_rna in source:
        gene_list.append(gene)
        ncpr_ps.append(logrank_ncpr)
        uversky_ps.append(logrank_uversky)
        
        ncpr_comp_rna_list.append(ncpr_comp_rna)
        ncpr_noncomp_rna_list.append(ncpr_noncomp_rna)
        uv_comp_rna_list.append(uv_comp_rna)
        uv_noncomp_rna_list.append(uv_noncomp_rna)

        
    
    results["gene"] = gene_list
    results["ncpr_p"] = ncpr_ps
    results["uversky_p"] = uversky_ps
    
    
    results["ncpr_comp_rna_exp"] = ncpr_comp_rna_list
    results["ncpr_noncomp_rna_rna_exp"] = ncpr_noncomp_rna_list
    results["ncpr_rna_comp_minus_noncomp"] = results["ncpr_comp_rna_exp"] - results["ncpr_noncomp_rna_rna_exp"]
    results["uversky_comp_rna_exp"] = uv_comp_rna_list
    results["uversky_noncomp_rna_rna_exp"] = uv_noncomp_rna_list
    results["uversky_rna_comp_minus_noncomp"] = results["uversky_comp_rna_exp"] - results["uversky_noncomp_rna_rna_exp"]
        
    return results

df = get_data()
results = rna_comp_corr()
results.to_csv("", index=False)  
    