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
        
        survivaltime = "OS_MONTHS"
        status = "OS_STATUS"
        
        logrank_ncpr = logrank_test(df1[survivaltime], df2[survivaltime], event_observed_A = df1[status], event_observed_B = df2[status])
        
        logrank_uversky = logrank_test(df3[survivaltime], df4[survivaltime], event_observed_A = df3[status], event_observed_B = df4[status])
        
        return gene, logrank_ncpr.p_value, logrank_uversky.p_value
        
    except:
        print 'error'
        result = pd.DataFrame()
        return result

genes = []
ncpr_ps = []
uversky_ps = []

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
    for gene, logrank_ncpr, logrank_uversky in source:
        genes.append(gene)
        ncpr_ps.append(logrank_ncpr)
        uversky_ps.append(logrank_uversky)
    
    results["gene"] = genes
    results["ncpr_p"] = ncpr_ps
    results["uversky_p"] = uversky_ps

df = get_data()
results = rna_comp_corr()
results.to_csv("")  
    