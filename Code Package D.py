import pandas as pd
from scipy.stats.stats import pearsonr
from lifelines import CoxPHFitter
from multiprocessing import Pool, cpu_count

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
        u_corr, u_p_value = pearsonr(df["CSU"], df[gene])
        n_corr, n_p_value = pearsonr(df["CSN"], df[gene])
        
        cox_ccu_os = coxcalc(df, gene, "OS_MONTHS", "OS_STATUS")
        cox_ccu_os["coef_os"] = cox_ccu_os["coef"]
        cox_ccu_os["p_os"] = cox_ccu_os["p"]
        cox_ccu_os = cox_ccu_os[["coef_os", "p_os"]]
        cox_ccu_os.insert(loc = 0, column = "gene", value = gene)
        cox_ccu_os.insert(loc = 1, column = "survival", value = "os")
        cox_ccu_os.insert(loc = 2, column = "os_n", value = len(df[gene]))
        
        cox_ccu_dfs = coxcalc(df, gene, "DFS_MONTHS", "DFS_STATUS")
    
        cox_ccu_os["coef_dfs"] = cox_ccu_dfs["coef"]
        cox_ccu_os["p_dfs"] = cox_ccu_dfs["p"]
        
        cox_ccu_os["n_pearson_r"] = n_corr
        cox_ccu_os["n_pearson_p"] = n_p_value
        cox_ccu_os["u_pearson_r"] = u_corr
        cox_ccu_os["u_pearson_p"] = u_p_value
        
        return cox_ccu_os
        
    except:
        print 'error'
        result = pd.DataFrame()
        return result

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
    for result in source:
        results = pd.concat([results, result], ignore_index=True)
    return results

df = get_data()
results = rna_comp_corr()
results.to_csv("")  
    