import pandas as pd
import numpy as np
from openpyxl import load_workbook
import os
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
from Bio import SeqIO
from localcider.sequenceParameters import SequenceParameters
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from multiprocessing import Pool, cpu_count
from sklearn import preprocessing


def aa_analysis(df, property):
    if property == "ncpr":
        df = df[pd.notnull(df['Amino_acids'])]
        df[["AA1","AA2"]] = df['Amino_acids'].str.split('/',expand=True)
        isoelectric_point = []
        for sequence in df["AA1"]:
            try:
                cdr3 = ProteinAnalysis(str(sequence))
                cidercdr3 = SequenceParameters(str(sequence)) 
                isoelectric_point.append(cidercdr3.get_NCPR())
            except:
                isoelectric_point.append(0)
                pass
        df["AA1_Iso"] = isoelectric_point
        isoelectric_point2 = []
        for sequence in df["AA2"]:
            try:
                cdr3 = ProteinAnalysis(str(sequence))
                cidercdr3 = SequenceParameters(str(sequence)) 
                isoelectric_point2.append(cidercdr3.get_NCPR())
            except:
                isoelectric_point2.append(0)
                pass
        df["AA2_Iso"] = isoelectric_point2
        df["AA_Iso_Delta"] = df["AA2_Iso"] - df["AA1_Iso"]
        df = df[["AA1_Iso", "AA2_Iso", "AA_Iso_Delta"]]
    elif property == "uversky_hydropathy":
        df = df[pd.notnull(df['Amino_acids'])]
        df[["AA1","AA2"]] = df['Amino_acids'].str.split('/',expand=True)
        isoelectric_point = []
        for sequence in df["AA1"]:
            try:
                cdr3 = ProteinAnalysis(str(sequence))
                cidercdr3 = SequenceParameters(str(sequence)) 
                isoelectric_point.append(cidercdr3.get_uversky_hydropathy())
            except:
                isoelectric_point.append(0)
                pass
        df["AA1_Iso"] = isoelectric_point
        isoelectric_point2 = []
        for sequence in df["AA2"]:
            try:
                cdr3 = ProteinAnalysis(str(sequence))
                cidercdr3 = SequenceParameters(str(sequence)) 
                isoelectric_point2.append(cidercdr3.get_uversky_hydropathy())
            except:
                isoelectric_point2.append(0)
                pass
        df["AA2_Iso"] = isoelectric_point2
        df["AA_Iso_Delta"] = df["AA2_Iso"] - df["AA1_Iso"]
        df = df[["AA1_Iso", "AA2_Iso", "AA_Iso_Delta"]]
    return df

def plotter(time_A, time_B, event_observed_A, event_observed_B):
    global property
    kmf = KaplanMeierFitter()
    kmf.fit(time_A, event_observed=event_observed_A, label="Top_%s"%property)
    ax = kmf.survival_function_.plot()
    kmf.fit(time_B, event_observed=event_observed_B, label="Bot_%s"%property)
    kmf.survival_function_.plot(ax=ax)
    fig = ax.get_figure()
    return fig

    
    
def getdata(cancer, receptor, sample, maindir):
    mutations = pd.read_csv(maindir + "VDJ Recoveries/{0}_Results/{1}.csv".format(cancer,"mutect"))
    
    
    mutations = mutations[mutations["Variant_Classification"] == "Missense_Mutation"]
    
    
    mutations["Tumor_Sample_Barcode"] = mutations["Tumor_Sample_Barcode"].apply(lambda x: x[0:12])
    cdr3 = pd.read_hdf(maindir+"VDJ Recoveries/{0}_Results/vdj_recoveries.h5".format(cancer), "physicochem")
    
    cdr3 = cdr3.loc[cdr3["Receptor"].str.contains(receptor)]
    
    cdr3 = cdr3.loc[cdr3["Sample"].str.contains(sample)]
    cdr3 = cdr3[["Filename", "ncpr", "uversky_hydropathy"]]
    cdr3["ncpr"] = cdr3["ncpr"].apply(pd.to_numeric, errors='coerce')
    cdr3["uversky_hydropathy"] = cdr3["uversky_hydropathy"].apply(pd.to_numeric, errors='coerce')
    cdr3max = cdr3.groupby("Filename", axis=0).max()
    cdr3min = cdr3.groupby("Filename", axis=0).min()
    cdr3 = cdr3.groupby("Filename", axis=0).mean()
    cdr3["ncpr_min"] = cdr3min["ncpr"]
    cdr3["ncpr_max"] = cdr3max["ncpr"]
    cdr3["uversky_hydropathy_min"] = cdr3min["uversky_hydropathy"]
    cdr3["uversky_hydropathy_max"] = cdr3max["uversky_hydropathy"]
    
    clinical = pd.read_csv(maindir+"VDJ Recoveries/%s_Results/clinical.csv"%cancer)
    clinical["OS_STATUS"] = np.where(clinical["OS_STATUS"]=="LIVING", 0, 1)
    clinical["DFS_STATUS"] = np.where(clinical["DFS_STATUS"]=="DiseaseFree", 0, 1)
    clinical = clinical.set_index("PATIENT_ID")

    mutations = mutations.join(cdr3, on='Tumor_Sample_Barcode')
    mutations = mutations[pd.notnull(mutations["ncpr"])]
    mutations = mutations[pd.notnull(mutations["ncpr_max"])]
    mutations = mutations[pd.notnull(mutations["ncpr_min"])]
    mutations = mutations[pd.notnull(mutations["uversky_hydropathy"])]
    mutations = mutations[pd.notnull(mutations["uversky_hydropathy_max"])]
    mutations = mutations[pd.notnull(mutations["uversky_hydropathy_min"])]
    mutations = mutations.join(clinical, on='Tumor_Sample_Barcode')
    mutations["Filename"] = mutations["Tumor_Sample_Barcode"]
    mutations = mutations.reset_index(drop=True)
    mutations[["ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]] = mutations[["ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]].apply(pd.to_numeric, errors='coerce')
    
    aa_analyzed = aa_analysis(mutations, "ncpr")
    mutations = mutations.join(aa_analyzed)
    aa_analyzed2 = aa_analysis(mutations, "uversky_hydropathy")
    aa_analyzed2["AA_Iso_Hydro"] = aa_analyzed2["AA_Iso_Delta"]
    mutations = mutations.join(aa_analyzed2["AA_Iso_Hydro"])
    mutants = mutations[['Filename', "Hugo_Symbol", "ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS", "AA_Iso_Delta", "AA_Iso_Hydro"]]
    
    cdr3 = cdr3.join(clinical)
    cdr3 = cdr3.reset_index()
    cdr3[["ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]] = cdr3[["ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]].apply(pd.to_numeric, errors='coerce')
    cdr3 = cdr3[['Filename', "ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]]
    return mutants, cdr3
   

def coxcalc(df, x, survivaltime, status):
    df5 = df[[status, survivaltime, x]]
    df5 = df5.dropna()
    cph = CoxPHFitter()
    cph.fit(df5, duration_col=survivaltime, event_col=status, show_progress=False)
    return cph.summary

def statscalc(x):
    print x
    global mutants
    global cdr3

    try:
        localmutants = mutants.loc[mutants["Hugo_Symbol"] == x]
        aalist = localmutants.groupby("Filename", axis=0).min()
        if len(aalist) >= 30:
            
            aalist["AA_Iso_Delta_Min"] = aalist["AA_Iso_Delta"]
            aalist["AA_Iso_Hydro_Min"] = aalist["AA_Iso_Hydro"]
            
            
            
            aalist2 = localmutants.groupby("Filename", axis=0).max()
            
            aalist2["AA_Iso_Delta_Max"] = aalist2["AA_Iso_Delta"]
            aalist2["AA_Iso_Hydro_Max"] = aalist2["AA_Iso_Hydro"]
            
            aalist3 = localmutants.groupby("Filename", axis=0).mean()
            
            localmutants = localmutants.drop_duplicates(subset="Filename")
            
            localmutants = localmutants.drop("AA_Iso_Delta", axis=1)
            localmutants = localmutants.drop("AA_Iso_Hydro", axis=1)
            
            localmutants = localmutants.join(aalist[["AA_Iso_Delta_Min", "AA_Iso_Hydro_Min"]], on="Filename")
    
            
            localmutants = localmutants.join(aalist2[["AA_Iso_Delta_Max", "AA_Iso_Hydro_Max"]], on="Filename")
            localmutants = localmutants.join(aalist3[["AA_Iso_Delta", "AA_Iso_Hydro"]], on="Filename")
            localcdr3 = cdr3.loc[~cdr3["Filename"].isin(localmutants["Filename"])]
            mutants_os = localmutants[["Filename", "Hugo_Symbol", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS", 'ncpr', 'ncpr_min', 'ncpr_max', 'uversky_hydropathy', 'uversky_hydropathy_min', 'uversky_hydropathy_max', "AA_Iso_Delta_Min", "AA_Iso_Delta_Max", "AA_Iso_Delta", "AA_Iso_Hydro_Min", "AA_Iso_Hydro_Max", "AA_Iso_Hydro"]]
            
            mutants_os = mutants_os.dropna(subset=["OS_STATUS", "OS_MONTHS", 'ncpr', 'ncpr_min', 'ncpr_max', 'uversky_hydropathy', 'uversky_hydropathy_min', 'uversky_hydropathy_max', "AA_Iso_Delta_Min", "AA_Iso_Delta_Max", "AA_Iso_Delta", "AA_Iso_Hydro_Min", "AA_Iso_Hydro_Max", "AA_Iso_Hydro"])
            mutants_os = mutants_os.reset_index(drop=True)
            mutants_dfs = localmutants[["DFS_STATUS", "DFS_MONTHS", 'ncpr', 'ncpr_min', 'ncpr_max', 'uversky_hydropathy', 'uversky_hydropathy_min', 'uversky_hydropathy_max', "AA_Iso_Delta_Min", "AA_Iso_Delta_Max", "AA_Iso_Delta", "AA_Iso_Hydro_Min", "AA_Iso_Hydro_Max", "AA_Iso_Hydro"]]
            mutants_dfs = mutants_dfs.dropna(subset=["DFS_STATUS", "DFS_MONTHS", 'ncpr', 'ncpr_min', 'ncpr_max', 'uversky_hydropathy', 'uversky_hydropathy_min', 'uversky_hydropathy_max', "AA_Iso_Delta_Min", "AA_Iso_Delta_Max", "AA_Iso_Delta", "AA_Iso_Hydro_Min", "AA_Iso_Hydro_Max", "AA_Iso_Hydro"])
            mutants_dfs = mutants_dfs.reset_index(drop=True)
        
        
            mutants_os["CS1N"] = mutants_os["ncpr_min"] * mutants_os["AA_Iso_Delta_Max"]
            mutants_os["CS2N"] = mutants_os["ncpr_max"] * mutants_os["AA_Iso_Delta_Min"]
            mutants_os["CSN"] =  mutants_os[["CS1N", "CS2N"]].min(axis=1)
            mutants_os["CSN"] = mutants_os["CSN"] * -1
            
    
            mutants_os["CSU"] = mutants_os["uversky_hydropathy_max"] * mutants_os["AA_Iso_Hydro_Max"]
            
            std_scale = preprocessing.StandardScaler().fit(mutants_os["CSN"].values.reshape(-1,1))
            mutants_os["CSN_scaled"] = std_scale.transform(mutants_os["CSN"].values.reshape(-1,1))
            
            std_scale = preprocessing.StandardScaler().fit(mutants_os["CSU"].values.reshape(-1,1))
            mutants_os["CSU_scaled"] = std_scale.transform(mutants_os["CSU"].values.reshape(-1,1))
            
            mutants_os["CS_sum"] = mutants_os["CSU_scaled"] + mutants_os["CSN_scaled"]
            mutants_os["CS_avg"] = mutants_os["CS_sum"]/2
            
    
            cox_ccu_os = coxcalc(mutants_os, "CSU", "OS_MONTHS", "OS_STATUS")
            cox_ccu_os.insert(loc = 0, column = "mutant", value = x)
            cox_ccu_os.insert(loc = 1, column = "tested", value = "uversky_hydropathy")
            cox_ccu_os.insert(loc = 2, column = "survival", value = "os")
            cox_ccu_os.insert(loc = 3, column = "n", value = len(mutants_os))
    
            cox_ccn_os = coxcalc(mutants_os, "CSN", "OS_MONTHS", "OS_STATUS")
            cox_ccn_os.insert(loc = 0, column = "mutant", value = x)
            cox_ccn_os.insert(loc = 1, column = "tested", value = "ncpr")
            cox_ccn_os.insert(loc = 1, column = "survival", value = "os")
            cox_ccn_os.insert(loc = 2, column = "n", value = len(mutants_os))
            
            cox_csa_os = coxcalc(mutants_os, "CS_avg", "OS_MONTHS", "OS_STATUS")
            cox_csa_os.insert(loc = 0, column = "mutant", value = x)
            cox_csa_os.insert(loc = 1, column = "tested", value = "CS_avg")
            cox_csa_os.insert(loc = 2, column = "survival", value = "os")
            cox_csa_os.insert(loc = 3, column = "n", value = len(mutants_os))
    
            
            mutants_dfs["CS1N"] = mutants_dfs["ncpr_min"] * mutants_dfs["AA_Iso_Delta_Max"]
            mutants_dfs["CS2N"] = mutants_dfs["ncpr_max"] * mutants_dfs["AA_Iso_Delta_Min"]
            mutants_dfs["CSN"] =  mutants_dfs[["CS1N", "CS2N"]].min(axis=1)
            mutants_dfs["CSN"] = mutants_dfs["CSN"] * -1
            
    
            mutants_dfs["CSU"] = mutants_dfs["uversky_hydropathy_max"] * mutants_dfs["AA_Iso_Hydro_Max"]
            
            
            std_scale = preprocessing.StandardScaler().fit(mutants_dfs["CSN"].values.reshape(-1,1))
            mutants_dfs["CSN_scaled"] = std_scale.transform(mutants_dfs["CSN"].values.reshape(-1,1))
            
            std_scale = preprocessing.StandardScaler().fit(mutants_dfs["CSU"].values.reshape(-1,1))
            mutants_dfs["CSU_scaled"] = std_scale.transform(mutants_dfs["CSU"].values.reshape(-1,1))
            
            mutants_dfs["CS_sum"] = mutants_dfs["CSU_scaled"] + mutants_dfs["CSN_scaled"]
            mutants_dfs["CS_avg"] = mutants_dfs["CS_sum"]/2
            
    
            cox_ccu_dfs = coxcalc(mutants_dfs, "CSU", "DFS_MONTHS", "DFS_STATUS")
            cox_ccu_dfs.insert(loc = 0, column = "mutant", value = x)
            cox_ccu_dfs.insert(loc = 1, column = "tested", value = "uversky_hydropathy")
            cox_ccu_dfs.insert(loc = 2, column = "survival", value = "dfs")
            cox_ccu_dfs.insert(loc = 3, column = "n", value = len(mutants_dfs))
    
            cox_ccn_dfs = coxcalc(mutants_dfs, "CSN", "DFS_MONTHS", "DFS_STATUS")
            cox_ccn_dfs.insert(loc = 0, column = "mutant", value = x)
            cox_ccn_dfs.insert(loc = 1, column = "tested", value = "ncpr")
            cox_ccn_dfs.insert(loc = 1, column = "survival", value = "dfs")
            cox_ccn_dfs.insert(loc = 2, column = "n", value = len(mutants_dfs))
            
            cox_csa_dfs = coxcalc(mutants_dfs, "CS_avg", "DFS_MONTHS", "DFS_STATUS")
            cox_csa_dfs.insert(loc = 0, column = "mutant", value = x)
            cox_csa_dfs.insert(loc = 1, column = "tested", value = "CS_avg")
            cox_csa_dfs.insert(loc = 2, column = "survival", value = "dfs")
            cox_csa_dfs.insert(loc = 3, column = "n", value = len(mutants_dfs))
    
    
            results = pd.concat([cox_ccu_os, cox_ccn_os, cox_csa_os, cox_ccu_dfs, cox_ccn_dfs, cox_csa_dfs], ignore_index=True)
            
            
            ##CSN KM using 0 as cut off
    
            mutants_os_top1 = mutants_os[mutants_os["CSN"] > 0]
            mutants_os_bot1 = mutants_os[mutants_os["CSN"] <= 0]
            
            mutants_dfs_top1 = mutants_dfs[mutants_dfs["CSN"] > 0]
            mutants_dfs_bot1 = mutants_dfs[mutants_dfs["CSN"] <= 0]
            
            logrank_mutants_os1 = logrank_test(mutants_os_top1["OS_MONTHS"], mutants_os_bot1["OS_MONTHS"], event_observed_A = mutants_os_top1["OS_STATUS"], event_observed_B = mutants_os_bot1["OS_STATUS"])
                        
            logrank_mutants_dfs1 = logrank_test(mutants_dfs_top1["DFS_MONTHS"], mutants_dfs_bot1["DFS_MONTHS"], event_observed_A = mutants_dfs_top1["DFS_STATUS"], event_observed_B = mutants_dfs_bot1["DFS_STATUS"])
            
            results["CSN_km_0_p_os"] = logrank_mutants_os1.p_value
            results["CSN_km_0_p_dfs"] = logrank_mutants_dfs1.p_value
            results["CSN_km_n_above_0"] = len(mutants_os_top1)
            results["CSN_km_n_lessorequal_0"] = len(mutants_os_bot1)
            
            ##CSN KM using median as cut off
            
    
            mutants_os_top2 = mutants_os[mutants_os["CSN"] > mutants_os["CSN"].median()]
            mutants_os_bot2 = mutants_os[mutants_os["CSN"] <= mutants_os["CSN"].median()]
            
            mutants_dfs_top2 = mutants_dfs[mutants_dfs["CSN"] > mutants_os["CSN"].median()]
            mutants_dfs_bot2 = mutants_dfs[mutants_dfs["CSN"] <= mutants_os["CSN"].median()]
            
            logrank_mutants_os2 = logrank_test(mutants_os_top2["OS_MONTHS"], mutants_os_bot2["OS_MONTHS"], event_observed_A = mutants_os_top2["OS_STATUS"], event_observed_B = mutants_os_bot2["OS_STATUS"])
                        
            logrank_mutants_dfs2 = logrank_test(mutants_dfs_top2["DFS_MONTHS"], mutants_dfs_bot2["DFS_MONTHS"], event_observed_A = mutants_dfs_top2["DFS_STATUS"], event_observed_B = mutants_dfs_bot2["DFS_STATUS"])
            
            results["CSN_km_50percentile_p_os"] = logrank_mutants_os2.p_value
            results["CSN_km_50percentile_p_dfs"] = logrank_mutants_dfs2.p_value
            results["CSN_km_50_top"] = len(mutants_os_top2)
            results["CSN_km_50_bot"] = len(mutants_os_bot2)
            
            
            ##CSU KM using median as cut off
            
            mutants_os_top3 = mutants_os[mutants_os["CSU"] > mutants_os["CSU"].median()]
            mutants_os_bot3 = mutants_os[mutants_os["CSU"] <= mutants_os["CSU"].median()]
            
            mutants_dfs_top3 = mutants_dfs[mutants_dfs["CSU"] > mutants_dfs["CSU"].median()]
            mutants_dfs_bot3 = mutants_dfs[mutants_dfs["CSU"] <= mutants_dfs["CSU"].median()]
            
            logrank_mutants_os3 = logrank_test(mutants_os_top3["OS_MONTHS"], mutants_os_bot3["OS_MONTHS"], event_observed_A = mutants_os_top3["OS_STATUS"], event_observed_B = mutants_os_bot3["OS_STATUS"])
                        
            logrank_mutants_dfs3 = logrank_test(mutants_dfs_top3["DFS_MONTHS"], mutants_dfs_bot3["DFS_MONTHS"], event_observed_A = mutants_dfs_top3["DFS_STATUS"], event_observed_B = mutants_dfs_bot3["DFS_STATUS"])
            
            results["CSU_km_50percentile_p_os"] = logrank_mutants_os3.p_value
            results["CSU_km_50percentile_p_dfs"] = logrank_mutants_dfs3.p_value
            results["CSU_km_50_top"] = len(mutants_os_top3)
            results["CSU_km_50_bot"] = len(mutants_os_bot3)
            
            ##CS_avg KM using median as cut off
            
            mutants_os_top4 = mutants_os[mutants_os["CS_avg"] > mutants_os["CS_avg"].median()]
            mutants_os_bot4 = mutants_os[mutants_os["CS_avg"] <= mutants_os["CS_avg"].median()]
            
            mutants_dfs_top4 = mutants_dfs[mutants_dfs["CS_avg"] > mutants_os["CS_avg"].median()]
            mutants_dfs_bot4 = mutants_dfs[mutants_dfs["CS_avg"] <= mutants_os["CS_avg"].median()]
            
            logrank_mutants_os4 = logrank_test(mutants_os_top4["OS_MONTHS"], mutants_os_bot4["OS_MONTHS"], event_observed_A = mutants_os_top4["OS_STATUS"], event_observed_B = mutants_os_bot4["OS_STATUS"])
                        
            logrank_mutants_dfs4 = logrank_test(mutants_dfs_top4["DFS_MONTHS"], mutants_dfs_bot4["DFS_MONTHS"], event_observed_A = mutants_dfs_top4["DFS_STATUS"], event_observed_B = mutants_dfs_bot4["DFS_STATUS"])
            
            results["CS_avg_km_50percentile_p_os"] = logrank_mutants_os4.p_value
            results["CS_avg_km_50percentile_p_dfs"] = logrank_mutants_dfs4.p_value
            results["CS_avg_km_50_top"] = len(mutants_os_top4)
            results["CS_avg_km_50_bot"] = len(mutants_os_bot4)
            
    
            results = results.rename(columns={'p': 'cox_p'})
            mutants_os = mutants_os.rename(index=str, columns={"ncpr":"ncpr_avg", "uversky_hydropathy":"uversky_hydropathy_avg"})
            
            return results, mutants_os
        
        else:
            results = pd.DataFrame()
            localmutants = pd.DataFrame()
            return results, localmutants
    except:
        results = pd.DataFrame()
        localmutants = pd.DataFrame()
        return results, localmutants

def statscalc_patientsum(mutants, cdr3):
    try:  
        localmutants = mutants.drop_duplicates(subset="Filename")
        
        localcdr3 = cdr3.loc[~cdr3["Filename"].isin(localmutants["Filename"])]
        mutants_os = localmutants[["Filename", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS", "ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "AA_Iso_Delta", "AA_Iso_Hydro"]]
        
        mutants_os = mutants_os.dropna(subset=["OS_STATUS", "OS_MONTHS", "ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min",  "AA_Iso_Delta", "AA_Iso_Hydro"])
        mutants_os = mutants_os.reset_index(drop=True)
        mutants_dfs = localmutants[["Filename", "DFS_STATUS", "DFS_MONTHS", "ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "AA_Iso_Delta", "AA_Iso_Hydro"]]
        mutants_dfs = mutants_dfs.dropna(subset=["DFS_STATUS", "DFS_MONTHS", "ncpr", "ncpr_max", "ncpr_min", "uversky_hydropathy", "uversky_hydropathy_max", "uversky_hydropathy_min", "AA_Iso_Delta", "AA_Iso_Hydro"])
        mutants_dfs = mutants_dfs.reset_index(drop=True)
    
        mutants_os["CS1N"] = mutants_os["ncpr_min"] * mutants_os["AA_Iso_Delta"]
        mutants_os["CS2N"] = mutants_os["ncpr_max"] * mutants_os["AA_Iso_Delta"]
        mutants_os["CSN"] =  mutants_os[["CS1N", "CS2N"]].min(axis=1)
        mutants_os["CSN"] = mutants_os["CSN"] * -1
    
    
        mutants_os["CSU"] = mutants_os["uversky_hydropathy_max"] * mutants_os["AA_Iso_Hydro"]
    
        mutants_os = mutants_os.groupby("Filename", axis=0).mean()
    
        std_scale = preprocessing.StandardScaler().fit(mutants_os["CSN"].values.reshape(-1,1))
        mutants_os["CSN_scaled"] = std_scale.transform(mutants_os["CSN"].values.reshape(-1,1))
        
        std_scale = preprocessing.StandardScaler().fit(mutants_os["CSU"].values.reshape(-1,1))
        mutants_os["CSU_scaled"] = std_scale.transform(mutants_os["CSU"].values.reshape(-1,1))
        
        mutants_os["CS_sum"] = mutants_os["CSU_scaled"] + mutants_os["CSN_scaled"]
        mutants_os["CS_avg"] = mutants_os["CS_sum"]/2
        
    
        cox_ccu_os = coxcalc(mutants_os, "CSU", "OS_MONTHS", "OS_STATUS")
        cox_ccu_os.insert(loc = 0, column = "tested", value = "uversky_hydropathy")
        cox_ccu_os.insert(loc = 1, column = "survival", value = "os")
        cox_ccu_os.insert(loc = 2, column = "n", value = len(mutants_os))
    
        cox_ccn_os = coxcalc(mutants_os, "CSN", "OS_MONTHS", "OS_STATUS")
        cox_ccn_os.insert(loc = 0, column = "tested", value = "ncpr")
        cox_ccn_os.insert(loc = 1, column = "survival", value = "os")
        cox_ccn_os.insert(loc = 2, column = "n", value = len(mutants_os))
        
        cox_csa_os = coxcalc(mutants_os, "CS_avg", "OS_MONTHS", "OS_STATUS")
        cox_csa_os.insert(loc = 0, column = "tested", value = "CS_avg")
        cox_csa_os.insert(loc = 1, column = "survival", value = "os")
        cox_csa_os.insert(loc = 2, column = "n", value = len(mutants_os))
    
        
        mutants_dfs["CS1N"] = mutants_dfs["ncpr_min"] * mutants_dfs["AA_Iso_Delta"]
        mutants_dfs["CS2N"] = mutants_dfs["ncpr_max"] * mutants_dfs["AA_Iso_Delta"]
        mutants_dfs["CSN"] =  mutants_dfs[["CS1N", "CS2N"]].min(axis=1)
        mutants_dfs["CSN"] = mutants_dfs["CSN"] * -1
        
    
        mutants_dfs["CSU"] = mutants_dfs["uversky_hydropathy_max"] * mutants_dfs["AA_Iso_Hydro"]
        
        
        std_scale = preprocessing.StandardScaler().fit(mutants_dfs["CSN"].values.reshape(-1,1))
        mutants_dfs["CSN_scaled"] = std_scale.transform(mutants_dfs["CSN"].values.reshape(-1,1))
        
        std_scale = preprocessing.StandardScaler().fit(mutants_dfs["CSU"].values.reshape(-1,1))
        mutants_dfs["CSU_scaled"] = std_scale.transform(mutants_dfs["CSU"].values.reshape(-1,1))
        
        mutants_dfs["CS_sum"] = mutants_dfs["CSU_scaled"] + mutants_dfs["CSN_scaled"]
        mutants_dfs["CS_avg"] = mutants_dfs["CS_sum"]/2
        
    
        cox_ccu_dfs = coxcalc(mutants_dfs, "CSU", "DFS_MONTHS", "DFS_STATUS")
    
        cox_ccu_dfs.insert(loc = 0, column = "tested", value = "uversky_hydropathy")
        cox_ccu_dfs.insert(loc = 1, column = "survival", value = "dfs")
        cox_ccu_dfs.insert(loc = 2, column = "n", value = len(mutants_dfs))
    
        cox_ccn_dfs = coxcalc(mutants_dfs, "CSN", "DFS_MONTHS", "DFS_STATUS")
    
        cox_ccn_dfs.insert(loc = 0, column = "tested", value = "ncpr")
        cox_ccn_dfs.insert(loc = 1, column = "survival", value = "dfs")
        cox_ccn_dfs.insert(loc = 2, column = "n", value = len(mutants_dfs))
        
        cox_csa_dfs = coxcalc(mutants_dfs, "CS_avg", "DFS_MONTHS", "DFS_STATUS")
    
        cox_csa_dfs.insert(loc = 0, column = "tested", value = "CS_avg")
        cox_csa_dfs.insert(loc = 1, column = "survival", value = "dfs")
        cox_csa_dfs.insert(loc = 2, column = "n", value = len(mutants_dfs))
    
    
        results = pd.concat([cox_ccu_os, cox_ccn_os, cox_csa_os, cox_ccu_dfs, cox_ccn_dfs, cox_csa_dfs], ignore_index=True)
        
        
        ##CSN KM using 0 as cut off
    
        mutants_os_top1 = mutants_os[mutants_os["CSN"] > 0]
        mutants_os_bot1 = mutants_os[mutants_os["CSN"] <= 0]
        
        mutants_dfs_top1 = mutants_dfs[mutants_dfs["CSN"] > 0]
        mutants_dfs_bot1 = mutants_dfs[mutants_dfs["CSN"] <= 0]
        
        logrank_mutants_os1 = logrank_test(mutants_os_top1["OS_MONTHS"], mutants_os_bot1["OS_MONTHS"], event_observed_A = mutants_os_top1["OS_STATUS"], event_observed_B = mutants_os_bot1["OS_STATUS"])
                    
        logrank_mutants_dfs1 = logrank_test(mutants_dfs_top1["DFS_MONTHS"], mutants_dfs_bot1["DFS_MONTHS"], event_observed_A = mutants_dfs_top1["DFS_STATUS"], event_observed_B = mutants_dfs_bot1["DFS_STATUS"])
        
        results["CSN_km_0_p_os"] = logrank_mutants_os1.p_value
        results["CSN_km_0_p_dfs"] = logrank_mutants_dfs1.p_value
        results["CSN_km_n_above_0"] = len(mutants_os_top1)
        results["CSN_km_n_lessorequal_0"] = len(mutants_os_bot1)
        
        ##CSN KM using median as cut off
        
    
        mutants_os_top2 = mutants_os[mutants_os["CSN"] > mutants_os["CSN"].median()]
        mutants_os_bot2 = mutants_os[mutants_os["CSN"] <= mutants_os["CSN"].median()]
        
        mutants_dfs_top2 = mutants_dfs[mutants_dfs["CSN"] > mutants_os["CSN"].median()]
        mutants_dfs_bot2 = mutants_dfs[mutants_dfs["CSN"] <= mutants_os["CSN"].median()]
        
        logrank_mutants_os2 = logrank_test(mutants_os_top2["OS_MONTHS"], mutants_os_bot2["OS_MONTHS"], event_observed_A = mutants_os_top2["OS_STATUS"], event_observed_B = mutants_os_bot2["OS_STATUS"])
                    
        logrank_mutants_dfs2 = logrank_test(mutants_dfs_top2["DFS_MONTHS"], mutants_dfs_bot2["DFS_MONTHS"], event_observed_A = mutants_dfs_top2["DFS_STATUS"], event_observed_B = mutants_dfs_bot2["DFS_STATUS"])
        
        results["CSN_km_50percentile_p_os"] = logrank_mutants_os2.p_value
        results["CSN_km_50percentile_p_dfs"] = logrank_mutants_dfs2.p_value
        results["CSN_km_50_top"] = len(mutants_os_top2)
        results["CSN_km_50_bot"] = len(mutants_os_bot2)
        
        
        ##CSU KM using median as cut off
        
        mutants_os_top3 = mutants_os[mutants_os["CSU"] > mutants_os["CSU"].median()]
        mutants_os_bot3 = mutants_os[mutants_os["CSU"] <= mutants_os["CSU"].median()]
        
        mutants_dfs_top3 = mutants_dfs[mutants_dfs["CSU"] > mutants_dfs["CSU"].median()]
        mutants_dfs_bot3 = mutants_dfs[mutants_dfs["CSU"] <= mutants_dfs["CSU"].median()]
        
        logrank_mutants_os3 = logrank_test(mutants_os_top3["OS_MONTHS"], mutants_os_bot3["OS_MONTHS"], event_observed_A = mutants_os_top3["OS_STATUS"], event_observed_B = mutants_os_bot3["OS_STATUS"])
                    
        logrank_mutants_dfs3 = logrank_test(mutants_dfs_top3["DFS_MONTHS"], mutants_dfs_bot3["DFS_MONTHS"], event_observed_A = mutants_dfs_top3["DFS_STATUS"], event_observed_B = mutants_dfs_bot3["DFS_STATUS"])
        
        results["CSU_km_50percentile_p_os"] = logrank_mutants_os3.p_value
        results["CSU_km_50percentile_p_dfs"] = logrank_mutants_dfs3.p_value
        results["CSU_km_50_top"] = len(mutants_os_top3)
        results["CSU_km_50_bot"] = len(mutants_os_bot3)
        
        ##CS_avg KM using median as cut off
        
        mutants_os_top4 = mutants_os[mutants_os["CS_avg"] > mutants_os["CS_avg"].median()]
        mutants_os_bot4 = mutants_os[mutants_os["CS_avg"] <= mutants_os["CS_avg"].median()]
        
        mutants_dfs_top4 = mutants_dfs[mutants_dfs["CS_avg"] > mutants_os["CS_avg"].median()]
        mutants_dfs_bot4 = mutants_dfs[mutants_dfs["CS_avg"] <= mutants_os["CS_avg"].median()]
        
        logrank_mutants_os4 = logrank_test(mutants_os_top4["OS_MONTHS"], mutants_os_bot4["OS_MONTHS"], event_observed_A = mutants_os_top4["OS_STATUS"], event_observed_B = mutants_os_bot4["OS_STATUS"])
                    
        logrank_mutants_dfs4 = logrank_test(mutants_dfs_top4["DFS_MONTHS"], mutants_dfs_bot4["DFS_MONTHS"], event_observed_A = mutants_dfs_top4["DFS_STATUS"], event_observed_B = mutants_dfs_bot4["DFS_STATUS"])
        
        results["CS_avg_km_50percentile_p_os"] = logrank_mutants_os4.p_value
        results["CS_avg_km_50percentile_p_dfs"] = logrank_mutants_dfs4.p_value
        results["CS_avg_km_50_top"] = len(mutants_os_top4)
        results["CS_avg_km_50_bot"] = len(mutants_os_bot4)
        
    
        results = results.rename(columns={'p': 'cox_p'})
        
        mutants_os = mutants_os.rename(index=str, columns={"ncpr":"ncpr_avg", "uversky_hydropathy":"uversky_hydropathy_avg"})
        
        return results, mutants_os

    except:
        results = pd.DataFrame()
        localmutants = pd.DataFrame()
        return results, localmutants
    
    
    
    

def kmcurve(mutants, cdr3):
    results = pd.DataFrame()
    raws = pd.DataFrame()
    p = Pool(cpu_count())
    source = p.map(statscalc, list(set(mutants["Hugo_Symbol"])))
    #for x in list(set(mutants["Hugo_Symbol"])): 
    for result, raw in source:
        #result, raw = statscalc(x) 
        results = pd.concat([results, result], ignore_index=True)
        raws = pd.concat([raws, raw], ignore_index=True)
    p.close()
    p.join()
    results = results.dropna()
    try:
        results = results.sort_values("cox_p", ascending=True)
    except:
        return results, raws
    return results, raws

maindir = ""
#maindir = ""
cancers = [ name for name in os.listdir(maindir + "VDJ Recoveries/") if os.path.isdir(os.path.join(maindir, "VDJ Recoveries/", name))]
print cancers
allresults= pd.DataFrame()
ptallresults= pd.DataFrame()
cols = ["cancer","sample","receptor","mutant","tested","survival","n","coef","exp(coef)","se(coef)","z","cox_p","lower 0.95","upper 0.95","CSN_km_0_p_os","CSN_km_0_p_dfs","CSN_km_n_above_0","CSN_km_n_lessorequal_0","CSN_km_50percentile_p_os","CSN_km_50percentile_p_dfs","CSN_km_50_top","CSN_km_50_bot","CSU_km_50percentile_p_os","CSU_km_50percentile_p_dfs","CSU_km_50_top","CSU_km_50_bot","CS_avg_km_50percentile_p_os","CS_avg_km_50percentile_p_dfs","CS_avg_km_50_top","CS_avg_km_50_bot"]
ptcols = ["cancer","sample","receptor","tested","survival","n","coef","exp(coef)","se(coef)","z","cox_p","lower 0.95","upper 0.95","CSN_km_0_p_os","CSN_km_0_p_dfs","CSN_km_n_above_0","CSN_km_n_lessorequal_0","CSN_km_50percentile_p_os","CSN_km_50percentile_p_dfs","CSN_km_50_top","CSN_km_50_bot","CSU_km_50percentile_p_os","CSU_km_50percentile_p_dfs","CSU_km_50_top","CSU_km_50_bot","CS_avg_km_50percentile_p_os","CS_avg_km_50percentile_p_dfs","CS_avg_km_50_top","CS_avg_km_50_bot"]
for cancer in cancers:
    for receptor in ["TRA|TRB", "TRA", "TRB", "IGH", "IGH|IGK|IGL", "IGH|IGK", "IGH|IGL"]:
        for sample in ["01|06|10", "01|06", "10"]:
            try:
                cancer = cancer.replace("_Results", "")
                mutationalgo = 'mutect' 
                ndiv = 2
                key = cancer+"_"+sample+"_"+receptor
                
                mutants, cdr3 = getdata(cancer, receptor, sample, maindir)
                
                
                ##Patient Scan
                
                ptresults, ptraws = statscalc_patientsum(mutants, cdr3)
                ptraws.to_hdf(maindir+"Complementarity/V3/missense/patient_rawV3.h5", key)
                ptraws.to_csv(maindir+"Complementarity/V3/missense/patient_raw_csvs/%s.csv"%key, index=True)
                
                ptresults.insert(0, "receptor", receptor)
                ptresults.insert(0, "sample", sample)
                ptresults.insert(0, "cancer", cancer)

                ptallresults = pd.concat([ptallresults, ptresults], ignore_index=True)
                
                ptresults = ptresults[ptcols]
                ptallresults = ptallresults[ptcols]
    

                
                ##Gene Scan
                
    
                results, raws = kmcurve(mutants, cdr3)
                raws.to_hdf(maindir+"Complementarity/V3/missense/gene_rawV3.h5", key)
                raws.to_csv(maindir+"Complementarity/V3/missense/gene_raw_csvs/%s.csv"%key, index=False)
                
                results.insert(0, "receptor", receptor)
                results.insert(0, "sample", sample)
                results.insert(0, "cancer", cancer)
    
                allresults = pd.concat([allresults, results], ignore_index=True)
    

                results = results[cols]
                allresults = allresults[cols]
    
                results.to_hdf(maindir+"Complementarity/V3/missense/gene_resultsV3.h5", key)

            
            except:
                print "ERROR ", cancer
                pass
            
def super_filter(allresults):
    
    filtered = allresults[allresults['cox_p'] < 0.05]
    filtered = filtered.sort_values("cox_p")
    
    super_filtered1 = filtered[filtered["tested"] == "ncpr"]
    super_filtered1 = super_filtered1[super_filtered1["survival"] == "os"]
    super_filtered1 = super_filtered1[super_filtered1["CSN_km_0_p_os"] <0.05]
    
    super_filtered2 = filtered[filtered["tested"] == "ncpr"]
    super_filtered2 = super_filtered2[super_filtered2["survival"] == "dfs"]
    super_filtered2 = super_filtered2[super_filtered2["CSN_km_0_p_dfs"] <0.05]
    
    super_filtered3 = filtered[filtered["tested"] == "uversky_hydropathy"]
    super_filtered3 = super_filtered3[super_filtered3["survival"] == "os"]
    super_filtered3 = super_filtered3[super_filtered3["CSU_km_50percentile_p_os"] <0.05]
    
    super_filtered4 = filtered[filtered["tested"] == "uversky_hydropathy"]
    super_filtered4 = super_filtered4[super_filtered4["survival"] == "dfs"]
    super_filtered4 = super_filtered4[super_filtered4["CSU_km_50percentile_p_dfs"] <0.05]
    
    super_filtered5 = filtered[filtered["tested"] == "CS_avg"]
    super_filtered5 = super_filtered5[super_filtered5["survival"] == "os"]
    super_filtered5 = super_filtered5[super_filtered5["CS_avg_km_50percentile_p_os"] <0.05]
    
    super_filtered6 = filtered[filtered["tested"] == "CS_avg"]
    super_filtered6 = super_filtered6[super_filtered6["survival"] == "dfs"]
    super_filtered6 = super_filtered6[super_filtered6["CS_avg_km_50percentile_p_dfs"] <0.05]
    
    super_filtered = pd.concat([super_filtered1, super_filtered2, super_filtered3, super_filtered4, super_filtered5, super_filtered6], ignore_index=True)
    
    return super_filtered

writer = pd.ExcelWriter(maindir+"Complementarity/V3/missense/gene_resultsV3.xlsx")
super_filtered = super_filter(allresults)
allresults = allresults[cols]
super_filtered = super_filtered[cols]
super_filtered.to_excel(writer, "filtered", index=False)
allresults.to_excel(writer, "all_results", index=False)
writer.save()

writer2 = pd.ExcelWriter(maindir+"Complementarity/V3/missense/patient_resultsV3.xlsx")
ptsuper_filtered = super_filter(ptallresults)
ptsuper_filtered.to_excel(writer2, "filtered", index=False)
ptallresults.to_excel(writer2, "all_results", index=False)
writer2.save()