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
from scipy import stats


class VDJrecord(object):
    
    def __init__(self, cancer, samples_path = None, vdjdb_path = None, clinical_path = None):
        self.cancer = cancer
        self.samples_path = samples_path
        self.vdjdb_path = vdjdb_path
        self.clinical_path = clinical_path
        self.projectid = "TCGA-%s" %self.cancer
    
    def load_raw_excel_path(self, path):
        print "Loading excel files from path: %s" %path
        self.raw = pd.DataFrame()
        for file in os.listdir(path):
            if file.endswith(".xlsx") and not file.startswith('~') and not file.startswith('V'):
                filepath = os.path.join(path, file)
                print "Loading: %s"%filepath
                df = pd.read_excel(filepath, sheetname=0)
                df["Receptor"] = file.replace(".xlsx", "")
                self.raw = pd.concat([self.raw, df])
        self.filtered = self.raw
        print "Loading excel files from path complete"
        return self.raw
        
    def load_raw_hdf_path(self, path):
        print "Loading hdf files from path: %s" %path
        self.raw = pd.DataFrame()
        for file in os.listdir(path):
            if file.endswith(".h5") and not file.startswith('~') and not file.startswith('v'):
                filepath = os.path.join(path, file)
                print "Loading: %s"%filepath
                df = pd.read_hdf(filepath, "raw")
                df["Receptor"] = file.replace(".h5", "")
                self.raw = pd.concat([self.raw, df])
        self.filtered = self.raw
        print "Loading hdf files from path complete"
        return self.raw
    
    def save_hdf(self, filepath):
        print "Saving data to hdf file: %s" %filepath
        if hasattr(self, 'raw'):
            print 'Saving raw data to hdf'
            self.raw.to_hdf(filepath, "raw")
        if hasattr(self, 'filtered_productive'):
            print 'Saving filtered productive data to hdf'
            self.filtered_productive.to_hdf(filepath, "filtered_productive")
        if hasattr(self, 'filtered_unproductive'):
            print 'Saving filtered unproductive data to hdf'
            self.filtered_unproductive.to_hdf(filepath, "filtered_unproductive")
        if hasattr(self, 'physicochem'):
            print 'Saving physicochemical data to hdf'
            self.physicochem.to_hdf(filepath, "physicochem")
        if hasattr(self, 'vdjdbmatch'):
            print 'Saving VDJDB matches to hdf'
            self.vdjdbmatch.to_hdf(filepath, "vdjdbmatch")
            
    def save_excel(self, filepath):
        print "Saving data to excel file: %s" %filepath
        writer = pd.ExcelWriter(filepath)
        if hasattr(self, 'raw'):
            print 'Saving raw data to excel'
            self.raw.to_excel(writer, "raw", index=False)
        if hasattr(self, 'filtered_productive'):
            print 'Saving filtered productive data to excel'
            self.filtered_productive.to_excel(writer, "filtered_productive", index=False)
        if hasattr(self, 'filtered_unproductive'):
            print 'Saving filtered unproductive data to excel'
            self.filtered_unproductive.to_excel(writer, "filtered_unproductive", index=False)
        if hasattr(self, 'physicochem'):
            print 'Saving physicochemical data to excel'
            self.physicochem.to_excel(writer, "physicochem", index=False)
        if hasattr(self, 'vdjdbmatch'):
            print 'Saving VDJDB matches to excel'
            self.vdjdbmatch.to_excel(writer, "vdjdbmatch", index=False)
        writer.save()

    def load_hdf(self, filepath):
        print "Loading data from hdf file: %s" %filepath
        try:
            print "Loading raw data from hdf file"
            self.raw = pd.read_hdf(filepath, "raw")
            self.filtered = self.raw
        except:
            raise ValueError("Hdf raw file doesn't exist")
        try:
            print "Loading filtered productive data from hdf file"
            self.filtered_productive = pd.read_hdf(filepath, "filtered_productive")
        except:
            print "Filtered productive data doesn't exist in file"
        try:
            print "Loading filtered unproductive data from hdf file"
            self.filtered_unproductive = pd.read_hdf(filepath, "filtered_unproductive")
        except:
            print "Filtered unproductive data doesn't exist in file"
        try:
            print "Loading physicochemical data from hdf file"
            self.physicochem = pd.read_hdf(filepath, "physicochem")
        except:
            print "Physicochemical data doesn't exist in file"
        try:
            print "Loading VDJDB matches from hdf file"
            self.vdjdbmatch = pd.read_hdf(filepath, "vdjdbmatch")
        except:
            print "VDJDB matches don't exist in file"


    def nt_match_length_filter(self, V_length = 19, J_length = 19):
        print "Filtering match lengths: V = {0}, J = {1}".format(V_length, J_length)
        self.filtered = self.filtered.loc[self.filtered["V Match Length"] >= V_length]
        self.filtered = self.filtered.loc[self.filtered["J Match Length"] >= J_length]
        return self.filtered

    def nt_match_percent_filter(self, V_percent = 90, J_percent = 90):
        print "Filtering match percents: V = {0}%, J = {1}%".format(V_percent, J_percent)
        self.filtered = self.filtered.loc[self.filtered["V Match Percent"] >= V_percent]
        self.filtered = self.filtered.loc[self.filtered["J Match Percent"] >= J_percent]
        return self.filtered
    
    def filename_format(self):
        if self.samples_path == None:
            raise ValueError("Must define path to samples db")
        else:
            print "Loading samples db: %s" %self.samples_path
            self.samples = pd.read_csv(self.samples_path)
            self.samples = self.samples[self.samples["Project ID"] == self.projectid]
            filenamelist = []
            samplelist = []
            sampletypelist = []
            self.filtered["Filename"] = self.filtered["Filename"].str.replace("sliced_","")
            self.filtered["Filename"] = self.filtered["Filename"].str.replace(".tsv","")
            print "Matching filenames to samples db"
            for i in self.filtered["Filename"]:
                filename = self.samples.loc[self.samples["File Name"] == i, ["Case ID", "Sample ID", "Sample Type"]]
                filenamelist.append(filename["Case ID"].iloc[0])
                samplelist.append(filename["Sample ID"].iloc[0][-3:])
                sampletypelist.append(filename["Sample Type"].iloc[0])
            self.filtered["Filename"] = filenamelist
            self.filtered.insert(1, 'Sample', samplelist)
            self.filtered.insert(2, 'Sample Type', sampletypelist)
            return self.filtered
    
    def raw_filename_format(self):
        if self.samples_path == None:
            raise ValueError("Must define path to samples db")
        else:
            print "Loading samples db: %s" %self.samples_path
            self.samples = pd.read_csv(self.samples_path)
            self.samples = self.samples[self.samples["Project ID"] == self.projectid]
            filenamelist = []
            samplelist = []
            sampletypelist = []
            self.raw["Filename"] = self.raw["Filename"].str.replace("sliced_","")
            self.raw["Filename"] = self.raw["Filename"].str.replace(".tsv","")
            print "Matching filenames to samples db"
            for i in self.raw["Filename"]:
                filename = self.samples.loc[self.samples["File Name"] == i, ["Case ID", "Sample ID", "Sample Type"]]
                filenamelist.append(filename["Case ID"].iloc[0])
                samplelist.append(filename["Sample ID"].iloc[0][-3:])
                sampletypelist.append(filename["Sample Type"].iloc[0])
            self.raw["Filename"] = filenamelist
            self.raw.insert(1, 'Sample', samplelist)
            self.raw.insert(2, 'Sample Type', sampletypelist)
            self.filtered = self.raw
            return self.raw
    
    def productive_unproductive_split(self):
        print "Splitting filtered data into productive and unproductive sets"
        self.filtered_unproductive = self.filtered[self.filtered["CDR3"] == 'Unproductive']
        self.filtered_productive = self.filtered[self.filtered["CDR3"] != 'Unproductive']
        return self.filtered_productive, self.filtered_unproductive
    
    def full_filter(self):
        self.raw_filename_format()
        self.nt_match_length_filter()
        self.nt_match_percent_filter()
        #self.filename_format()
        self.productive_unproductive_split()
        return self.filtered_productive, self.filtered_unproductive
    
    def vdjdb_match(self):
        if self.vdjdb_path == None:
            raise ValueError("Must define path to samples db")
        else:
            print "Matching productive CDR3s to VDJDB"
            tra = self.raw.loc[self.raw["Receptor"] == "TRA"]
            trb = self.raw.loc[self.raw["Receptor"] == "TRB"]
            vdjdb = pd.read_csv(self.vdjdb_path)
            travdjdb = vdjdb[vdjdb["Gene"] == "TRA"]
            trbvdjdb = vdjdb[vdjdb["Gene"] == "TRB"]
            tra = tra.merge(travdjdb, on='CDR3', how='inner')
            trb = trb.merge(trbvdjdb, on='CDR3', how='inner')
            df = pd.concat([tra, trb])
            df = df.drop_duplicates("Read ID")
            df = df.drop(["complex.id", "Gene", "V", "J"], axis=1)
            self.vdjdbmatch = df
            return self.vdjdbmatch
    
    def analyze_physicochem(self):
        print "Analyzing CDR3 physicochemical data"
        df = self.filtered_productive[["Filename", "CDR3", "Sample", "Sample Type", "Receptor"]]
        length = []
        fraction_tiny=[]
        fraction_small=[]
        fraction_charged=[]
        fraction_positive=[]
        fraction_negative=[]
        fraction_expanding=[]
        fraction_aromatic=[]
        molecular_weight = []
        isoelectric_point = []
        gravy = []
        aromaticity = []
        instability_index = []
        secondary_structure_helix = []
        secondary_structure_turn = []
        secondary_structure_sheet = []
        mean_hydropathy = []
        uversky_hydropathy = []
        NCPR = []
        kappa = []
        omega = []
        PPII_propensity = []
        delta = []
        fraction_disorder_promoting = []
        
        def count_aa(sequence, aas):
            count = 0
            for aa in aas:
                count += sequence.count(aa)
            return count
            
        def sequence_analysis(sequence):
            try:
                cdr3 = ProteinAnalysis(str(sequence))
                cidercdr3 = SequenceParameters(str(sequence)) 
                molecular_weight = (cdr3.molecular_weight())
                isoelectric_point = (cdr3.isoelectric_point())
                gravy = (cdr3.gravy())
                aromaticity = (cdr3.aromaticity())
                instability_index = (cdr3.instability_index())
                secondary = cdr3.secondary_structure_fraction()
                secondary_structure_helix = (secondary[0])
                secondary_structure_turn = (secondary[1])
                secondary_structure_sheet = (secondary[2])
                length = (len(sequence))
                fraction_tiny = (float(count_aa(sequence,"ABCGST"))/float(len(sequence)))
                fraction_small = (float(count_aa(sequence,"ABCDGNPSTV"))/float(len(sequence)))
                fraction_aromatic = (float(count_aa(sequence,"FHWY"))/float(len(sequence)))
                fraction_charged = (cidercdr3.get_FCR())
                fraction_positive = (cidercdr3.get_fraction_positive())
                fraction_negative = (cidercdr3.get_fraction_negative())
                fraction_expanding = (cidercdr3.get_fraction_expanding())
                mean_hydropathy = (cidercdr3.get_mean_hydropathy())
                NCPR = (cidercdr3.get_NCPR())
                uversky_hydropathy = (cidercdr3.get_uversky_hydropathy())
                kappa = (cidercdr3.get_mean_hydropathy())
                omega = (cidercdr3.get_Omega())
                PPII_propensity = (cidercdr3.get_PPII_propensity())
                delta = (cidercdr3.get_delta())
                fraction_disorder_promoting = (cidercdr3.get_fraction_disorder_promoting())
                return (molecular_weight, isoelectric_point, gravy, aromaticity, instability_index, secondary_structure_helix, secondary_structure_turn, secondary_structure_sheet, length, fraction_tiny, fraction_small, fraction_aromatic, fraction_charged, fraction_positive, fraction_negative, fraction_expanding, mean_hydropathy, NCPR, uversky_hydropathy,  kappa, omega, PPII_propensity, delta, fraction_disorder_promoting)
            except:
                return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,  np.nan, np.nan, np.nan, np.nan, np.nan)

        for sequence in df["CDR3"]:
            a = sequence_analysis(sequence)
            molecular_weight.append(a[0])
            isoelectric_point.append(a[1])
            gravy.append(a[2])
            aromaticity.append(a[3])
            instability_index.append(a[4])
            secondary_structure_helix.append(a[5])
            secondary_structure_turn.append(a[6])
            secondary_structure_sheet.append(a[7])
            length.append(a[8])
            fraction_tiny.append(a[9])
            fraction_small.append(a[10])
            fraction_aromatic.append(a[11])
            fraction_charged.append(a[12])
            fraction_positive.append(a[13])
            fraction_negative.append(a[14])
            fraction_expanding.append(a[15])
            mean_hydropathy.append(a[16])
            NCPR.append(a[17])
            uversky_hydropathy.append(a[18])
            kappa.append(a[19])
            omega.append(a[20])
            PPII_propensity.append(a[21])
            delta.append(a[22])
            fraction_disorder_promoting.append(a[23])
        df["length"] = length
        df["fraction_tiny"] = fraction_tiny
        df["fraction_small"] = fraction_small
        df["fraction_aromatic"] = fraction_aromatic
        df["fraction_charged"] = fraction_charged
        df["fraction_positive"] = fraction_positive
        df["fraction_negative"] = fraction_negative
        df["fraction_expanding"] = fraction_expanding
        df["fraction_disorder_promoting"] = fraction_disorder_promoting
        df["molecular_weight"] = molecular_weight
        df["isoelectric_point"] = isoelectric_point
        df["ncpr"] = NCPR
        df["mean_hydropathy"] = mean_hydropathy
        df["uversky_hydropathy"] = uversky_hydropathy
        df["ppii_propensity"] = PPII_propensity
        df["gravy"] = gravy
        df["aromaticity"] = aromaticity
        df["kappa"] = kappa
        df["omega"] = omega
        df["delta"] = delta
        df["instability_index"] = instability_index
        df["secondary_structure_helix"] = secondary_structure_helix
        df["secondary_structure_turn"] = secondary_structure_turn
        df["secondary_structure_sheet"] = secondary_structure_sheet
        self.physicochem = df
        return self.physicochem
    
    
class VDJanalyze(object):
    
    def __init__(self, vdjrecord, rootdir):
        self.productive = vdjrecord.filtered_productive
        self.physicochem = vdjrecord.physicochem
        self.vdjdbmatch = vdjrecord.vdjdbmatch
        self.rootdir = rootdir
        clinical_path = os.path.join(self.rootdir,"clinical.csv")
        print "Loading clinical data from: %s" %clinical_path
        self.clinical = pd.read_csv(clinical_path)
        self.clinical["Filename"] = self.clinical["PATIENT_ID"]
        self.clinical["OS_STATUS"] = np.where(self.clinical["OS_STATUS"]=="LIVING", 0, 1)
        self.clinical["DFS_STATUS"] = np.where(self.clinical["DFS_STATUS"]=="DiseaseFree", 0, 1)
        self.km_receptors_plots_directory = self.rootdir+"/km_receptors_plots/"
        self.km_physicochem_plots_directory = self.rootdir+"/km_physicochem_plots/"
        if not os.path.exists(self.km_receptors_plots_directory):
            os.makedirs(self.km_receptors_plots_directory)
        if not os.path.exists(self.km_physicochem_plots_directory):
            os.makedirs(self.km_physicochem_plots_directory)
            
    def save_hdf(self, filepath):
        print "Saving data to hdf file: %s" %filepath
        if hasattr(self, 'km_receptors_results'):
            print 'Saving receptor km data to hdf'
            self.km_receptors_results.to_hdf(filepath, "km_receptors")
        if hasattr(self, 'km_physicochem_results'):
            print 'Saving physicochemical km data to hdf'
            self.km_physicochem_results.to_hdf(filepath, "km_physicochem")
        if hasattr(self, 'cox_physicochem_results'):
            print 'Saving physicochemical cox regression data to hdf'
            self.cox_physicochem_results.to_hdf(filepath, "cox_physicochem")
        if hasattr(self, 'km_vdjdb_results'):
            print 'Saving vdjdb km data to hdf'
            self.km_vdjdb_results.to_hdf(filepath, "km_vdjdb")
        
            
    def save_excel(self, filepath):
        print "Saving data to excel file: %s" %filepath
        writer = pd.ExcelWriter(filepath)
        if hasattr(self, 'km_receptors_results'):
            print 'Saving receptor km data to excel'
            self.km_receptors_results.to_excel(writer, "km_receptors", index=False)
        if hasattr(self, 'km_physicochem_results'):
            print 'Saving physicochemical km data to excel'
            self.km_physicochem_results.to_excel(writer, "km_physicochem", index=False)
        if hasattr(self, 'cox_physicochem_results'):
            print 'Saving physicochemical cox regression data to excel'
            self.cox_physicochem_results.to_excel(writer, "cox_physicochem", index=False)
        if hasattr(self, 'km_vdjdb_results'):
            print 'Saving vdjdb km data to excel'
            self.km_vdjdb_results.to_excel(writer, "km_vdjdb", index=False)
        writer.save()
        
    def load_hdf(self, filepath):
        print "Loading data from hdf file: %s" %filepath
        try:
            print "Loading physicochem KM results from hdf file"
            self.km_physicochem_results = pd.read_hdf(filepath, "km_physicochem")
        except:
            print "Physicochem KM data doesn't exist in file"
        
        try:
            print "Loading physicochem cox regression results from hdf file"
            self.cox_physicochem_results = pd.read_hdf(filepath, "cox_physicochem")
        except:
            print "Physicochem cox regression data doesn't exist in file"
            
        try:
            print "Loading receptor KM results from hdf file"
            self.km_receptors_results = pd.read_hdf(filepath, "km_receptors")
        except:
            print "Receptors KM data doesn't exist in file"
            
        try:
            print "Loading vdjdb KM results from hdf file"
            self.km_vdjdb_results = pd.read_hdf(filepath, "km_vdjdb")
        except:
            print "vdjdb KM data doesn't exist in file"
        
    
    def plotter(self, time_A, time_B, event_observed_A, event_observed_B, label_A, label_B):
        kmf = KaplanMeierFitter()
        kmf.fit(time_A, event_observed=event_observed_A, label=label_A)
        ax = kmf.survival_function_.plot()
        kmf.fit(time_B, event_observed=event_observed_B, label=label_B)
        kmf.survival_function_.plot(ax=ax)
        ax.set_ylim([0,1])
        fig = ax.get_figure()
        return fig
    
    def medians(self, time_A, time_B, event_observed_A, event_observed_B, label_A, label_B):
        kmfa = KaplanMeierFitter()
        kmfb = KaplanMeierFitter()
        kmfa.fit(time_A, event_observed=event_observed_A)
        kmfb.fit(time_B, event_observed=event_observed_B)
        return kmfa.median_, kmfb.median_
        
    
    def km_vdjdb(self):
        def statscalc(epitope, receptor, site):
            try:
                localvdjdb = self.vdjdbmatch
                localvdjdb = localvdjdb[localvdjdb["Epitope species"] == epitope]
                localvdjdb = localvdjdb[localvdjdb["Receptor"].str.contains(receptor)]
                if site == "blood":
                    localvdjdb = localvdjdb[localvdjdb["Sample"].str.contains("10")]
                elif site == "tissue":
                    localvdjdb = localvdjdb[~localvdjdb["Sample"].str.contains("10")]
                
                localclinical = self.clinical[["Filename", "OS_STATUS", "OS_MONTHS", "DFS_MONTHS", "DFS_STATUS"]]
                localclinical = localclinical[localclinical.columns.difference(['Filename'])].apply(pd.to_numeric, errors='coerce')
                
                localclinical["Filename"] = self.clinical["Filename"]
                            
                localclinical_os = localclinical[["Filename", "OS_STATUS", "OS_MONTHS"]]
                localclinical_os = localclinical_os.dropna()
                
                localclinical_dfs = localclinical[["Filename", "DFS_STATUS", "DFS_MONTHS"]]
                
                localclinical_dfs = localclinical_dfs.dropna()
                
                
                df3 = localclinical_os[localclinical_os["Filename"].isin(localvdjdb["Filename"])]
                df4 = localclinical_os[~localclinical_os["Filename"].isin(localvdjdb["Filename"])]
                
                df5 = localclinical_dfs[localclinical_dfs["Filename"].isin(localvdjdb["Filename"])]
                df6 = localclinical_dfs[~localclinical_dfs["Filename"].isin(localvdjdb["Filename"])]
                    
    
                logrank_os = logrank_test(df3["OS_MONTHS"], df4["OS_MONTHS"], event_observed_A = df3["OS_STATUS"], event_observed_B = df4["OS_STATUS"])
                logrank_dfs = logrank_test(df5["DFS_MONTHS"], df6["DFS_MONTHS"], event_observed_A = df5["DFS_STATUS"], event_observed_B = df6["DFS_STATUS"])
    
                    
                median_res_os = self.medians(df3["OS_MONTHS"], df4["OS_MONTHS"], event_observed_A = df3["OS_STATUS"], event_observed_B = df4["OS_STATUS"], label_A = "Top", label_B = "Bottom")
                pos_survival_avg_os = float(median_res_os[0])
                neg_survival_avg_os  = float(median_res_os[1])
                
                
                
                median_res_dfs = self.medians(df5["DFS_MONTHS"], df6["DFS_MONTHS"], event_observed_A = df5["DFS_STATUS"], event_observed_B = df6["DFS_STATUS"], label_A = "Top", label_B = "Bottom")
                pos_survival_avg_dfs = float(median_res_dfs[0])
                neg_survival_avg_dfs  = float(median_res_dfs[1])    
                    
                return logrank_os.p_value, logrank_os.test_statistic, logrank_dfs.p_value, logrank_dfs.test_statistic, len(df3), len(df4), len(df5), len(df6), pos_survival_avg_os, neg_survival_avg_os, pos_survival_avg_dfs, neg_survival_avg_dfs
            
            except:
                return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        
        
        epitopelist = []
        receptorlist = []
        sitelist = []
        logrank_p_os = []
        logrank_f_os = []
        pos_survival_avg_os = []
        neg_survival_avg_os = []
        n_pos_os = []
        n_neg_os = []
        logrank_p_dfs = []
        logrank_f_dfs = []
        pos_survival_avg_dfs = []
        neg_survival_avg_dfs = []
        n_pos_dfs = []
        n_neg_dfs = []
        
        for receptor in ["TRA", "TRB", "TRA|TRB"]:
            for site in ["tissue", "blood", "tissue|blood"]:
                for epitope in list(self.vdjdbmatch["Epitope species"].unique()):
                    a = statscalc(epitope, receptor, site)
                    epitopelist.append(epitope)
                    receptorlist.append(receptor)
                    sitelist.append(site)
                    logrank_p_os.append(a[0])
                    logrank_f_os.append(a[1])
                    logrank_p_dfs.append(a[2])
                    logrank_f_dfs.append(a[3])
                    n_pos_os.append(a[4])
                    n_neg_os.append(a[5])
                    n_pos_dfs.append(a[6])
                    n_neg_dfs.append(a[7])
                    pos_survival_avg_os.append(a[8])
                    neg_survival_avg_os.append(a[9])
                    pos_survival_avg_dfs.append(a[10])
                    neg_survival_avg_dfs.append(a[11])
        
        results = pd.DataFrame()
        results["epitope"] = epitopelist
        results["site"] = sitelist
        results["receptor"] = receptorlist
        results["logrank_p_os"] = logrank_p_os
        results["logrank_f_os"] = logrank_f_os
        results["logrank_p_dfs"] = logrank_p_dfs
        results["logrank_f_dfs"] = logrank_f_dfs
        results["n_epitope_os"] = n_pos_os
        results["n_ar_os"] = n_neg_os
        results["n_epitope_dfs"] = n_pos_dfs
        results["n_ar_dfs"] = n_neg_dfs
        results["epitope_survival_med_os"] = pos_survival_avg_os
        results["ar_survival_med_os"] = neg_survival_avg_os
        results["epitope_survival_med_dfs"] = pos_survival_avg_dfs
        results["ar_survival_med_dfs"] = neg_survival_avg_dfs
        results = results.dropna()
        self.km_vdjdb_results = results.sort_values(by="logrank_p_os")
        self.km_vdjdb_results = self.km_vdjdb_results.reset_index(drop=True)
        return self.km_vdjdb_results

    def cox_physicochem(self):            
        def statscalc(x, df2, survivaltime, status):
            df5 = df2[[status, survivaltime, x]]
            cph = CoxPHFitter()
            cph.fit(df5, duration_col=survivaltime, event_col=status, show_progress=False)
            return cph.summary
        
        def kmcurve(df, df2, site, receptor, survivaltime, status):
            print "Analyzing CDR3 Physicochemical Cox Regression. Parameters: Receptor - {0}, Site - {1}, Survival - {2}".format(receptor, site, survivaltime.replace("_MONTHS",""))
            
            results = pd.DataFrame()
            df2 = df2[df2["Receptor"].str.contains(receptor)]
            for x in list(df.columns.values):
                result = statscalc(x, df2, survivaltime, status)
                results = pd.concat([results, result])
            results.insert(loc = 0, column = "site", value = site)
            results.insert(loc = 1, column = "receptor", value = receptor)
            results.insert(loc = 2, column = "survival_type", value = survivaltime.replace("_MONTHS", ""))
            results.insert(loc = 3, column = "n", value = len(df2))
            results = results.dropna()
            results = results.reset_index().rename(columns={results.index.name:'attribute'})
            return results

        localphysicochem = self.physicochem.set_index("Filename")
        localclinical = self.clinical[["Filename", "OS_STATUS", "OS_MONTHS", "DFS_MONTHS", "DFS_STATUS"]].join(localphysicochem, on='Filename')
        localclinical = localclinical[localclinical.columns.difference(["CDR3", "Sample Type"])]
        localclinical = localclinical.loc[:,~localclinical.columns.duplicated()]
        localphysicochem = localphysicochem.loc[:,~localphysicochem.columns.duplicated()]
        
        
        localclinical_dfs = localclinical[localclinical.columns.difference(['OS_MONTHS', 'OS_STATUS'])]
        localclinical = localclinical[localclinical.columns.difference(['DFS_MONTHS', 'DFS_STATUS'])]
        
        localclinical[localclinical.columns.difference(['Filename', 'Sample', "Receptor"])] = localclinical[localclinical.columns.difference(['Filename', 'Sample', "Receptor"])].apply(pd.to_numeric, errors='coerce')
        localclinical = localclinical.dropna()
        localclinical = localclinical.set_index(["Receptor", "Filename"])
        
        localclinical_dfs[localclinical_dfs.columns.difference(['Filename', 'Sample', "Receptor"])] = localclinical_dfs[localclinical_dfs.columns.difference(['Filename', 'Sample', "Receptor"])].apply(pd.to_numeric, errors='coerce')
        localclinical_dfs = localclinical_dfs.dropna()
        localclinical_dfs = localclinical_dfs.set_index(["Receptor", "Filename"])
        
        tissue = localclinical.loc[localclinical["Sample"].str.contains("10")==False]
        tissue = tissue.loc[tissue["Sample"].str.contains("11")==False]
        tissue = tissue.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        tissue = tissue.reset_index()
        blood = localclinical.loc[localclinical["Sample"].str.contains("10")]
        blood = blood.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        blood = blood.reset_index()
        
        
        tissue_dfs = localclinical_dfs.loc[localclinical_dfs["Sample"].str.contains("10")==False]
        tissue_dfs = tissue_dfs.loc[tissue_dfs["Sample"].str.contains("11")==False]
        tissue_dfs = tissue_dfs.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        tissue_dfs = tissue_dfs.reset_index()
        blood_dfs = localclinical_dfs.loc[localclinical_dfs["Sample"].str.contains("10")]
        blood_dfs = blood_dfs.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        blood_dfs = blood_dfs.reset_index()
        
        localphysicochem = localphysicochem.drop(["CDR3", "Sample", "Sample Type", "Receptor"], axis=1)
        self.cox_physicochem_results = pd.DataFrame()
        
        
        for receptor in ["TR", "TRA|TRB", "TRG|TRD", "IGH|IGL", "IGH|IGK", "IG", "TRA", "TRB", "TRD", "TRG", "IGK", "IGL", "IGH"]:
            for survivaltype in ["OS", "DFS"]:
                for site in ["TISSUE", "BLOOD"]:
                    if survivaltype == "OS":
                        if site == "TISSUE":
                            try:
                                results = kmcurve(localphysicochem, tissue, site, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype)
                                self.cox_physicochem_results = pd.concat([self.cox_physicochem_results, results])
                            except:
                                pass
                        if site == "BLOOD":
                            try:
                                results = kmcurve(localphysicochem, blood, site, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype)
                                self.cox_physicochem_results = pd.concat([self.cox_physicochem_results, results])
                            except:
                                pass
                                
                    if survivaltype == "DFS":
                        if site == "TISSUE":
                            try:
                                results = kmcurve(localphysicochem, tissue_dfs, site, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype)
                                self.cox_physicochem_results = pd.concat([self.cox_physicochem_results, results])
                            except:
                                pass

                        if site == "BLOOD":
                            try:
                                results = kmcurve(localphysicochem, blood_dfs, site, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype)
                                self.cox_physicochem_results = pd.concat([self.cox_physicochem_results, results])
                            except:
                                pass

        self.cox_physicochem_results = self.cox_physicochem_results.sort_values(by="p")
        self.cox_physicochem_results = self.cox_physicochem_results.reset_index(drop=True)
        return self.cox_physicochem_results
        


    def km_physicochem(self, plot=False):
        def statscalc(x, df2, receptor, site, ndiv, survivaltime, status, plot):
            df5 = df2[[status, survivaltime, x]]
            df5 = df5.dropna()
            df5 = df5.sort_values(by=[x])
            df5 = df5.reset_index()
            df4 = df5[:(len(df5)/ndiv)]
            df3 = df5[(len(df5)-(len(df5)/ndiv)):]
            logrank = logrank_test(df3[survivaltime], df4[survivaltime], event_observed_A = df3[status], event_observed_B = df4[status])
            if plot == True:
                fig = self.plotter(df3[survivaltime], df4[survivaltime], event_observed_A = df3[status], event_observed_B = df4[status], label_A = "Top %s"%x, label_B = "Bottom %s"%x)
                fig.savefig(self.km_physicochem_plots_directory+"/{0}_{1}_{2}_{3}_{4}.pdf".format(x, receptor, "DIV%s"%ndiv, survivaltime[:survivaltime.find("_")], site))
            n = len(df5)
            median_res = self.medians(df3[survivaltime], df4[survivaltime], event_observed_A = df3[status], event_observed_B = df4[status], label_A = "Top %s"%x, label_B = "Bottom %s"%x)
            top_survival_avg = float(median_res[0])
            bot_survival_avg  = float(median_res[1])
            ttest = stats.ttest_ind(df3[x],df4[x], equal_var=False)
            tf = ttest[0]
            tp = ttest[1]
            df3["logexpression"] = df3[x] + 1
            df4["logexpression"] = df4[x] + 1
            logtop = np.log10(df3["logexpression"].astype(float))
            logbot = np.log10(df4["logexpression"].astype(float))
            logt = stats.ttest_ind(logtop,logbot, equal_var=False)
            logtf = logt[0]
            logtp = logt[1]
            try:
                mwu = stats.mannwhitneyu(df3[x],df4[x])
            except:
                mwu = ('mwu_error', 'mwu_error')
            mwuf = mwu[0]
            mwup = mwu[1]
            top_expression_avg = df3[x].mean()
            bot_expression_avg = df4[x].mean()
            return (x, logrank.p_value, logrank.test_statistic, n, top_survival_avg, bot_survival_avg, tf, tp, logtf, logtp, mwuf, mwup, top_expression_avg, bot_expression_avg)
            
        def kmcurve(df, df2, site, ndiv, receptor, survivaltime, status, plot):
            print "Analyzing CDR3 Physicochemical KM Curves. Parameters: Receptor - {0}, Site - {1}, Ndiv - {2}, Survival - {3}".format(receptor, site, ndiv, survivaltime.replace("_MONTHS",""))
            genes = []
            logrank_f = []
            logrank_p = []
            n = []
            top_survival_avg = []
            bot_survival_avg = []
            tf = []
            tp = []
            logtf = []
            logtp = []
            mwuf = []
            mwup = []
            top_expression_avg = []
            bot_expression_avg = []
            df2 = df2[df2["Receptor"].str.contains(receptor)]
            for x in list(df.columns.values):
                a = statscalc(x, df2, receptor, site, ndiv, survivaltime, status, plot)
                genes.append(a[0])
                logrank_p.append(a[1])
                logrank_f.append(a[2])
                n.append(a[3])
                top_survival_avg.append(a[4])
                bot_survival_avg.append(a[5])
                tf.append(a[6])
                tp.append(a[7])
                logtf.append(a[8])
                logtp.append(a[9])
                mwuf.append(a[10])
                mwup.append(a[11])
                top_expression_avg.append(a[12])
                bot_expression_avg.append(a[13])
            results = pd.DataFrame()
            results["attribute"] = genes
            results["site"] = site
            results["receptor"] = receptor
            results["ndiv"] = ndiv
            results["survival_type"] = survivaltime.replace("_MONTHS", "")
            results["logrank_fvalue"] = logrank_f
            results["logrank_pvalue"] = logrank_p
            results["npatients"] = n
            results["top_survival_median"] = top_survival_avg
            results["bot_survival_median"] = bot_survival_avg
            results["survival_difference"] = results["top_survival_median"] - results["bot_survival_median"]
            results["expression_ttest_f"] = tf
            results["expression_ttest_p"] = tp
            results["expression_log_ttest_f"] = logtf
            results["expression_log_ttest_p"] = logtp
            results["expression_mwu_f"] = mwuf
            results["expression_mwu_p"] = mwup
            results["top_expression_avg"] = top_expression_avg
            results["bot_expression_avg"] = bot_expression_avg
            results = results.dropna()
            results = results.sort_values("logrank_pvalue", ascending=True)
            results = results.reset_index(drop = True)
            return results
            
        localphysicochem = self.physicochem.set_index("Filename")
        localclinical = self.clinical[["Filename", "OS_STATUS", "OS_MONTHS", "DFS_MONTHS", "DFS_STATUS"]].join(localphysicochem, on='Filename')
        localclinical = localclinical[localclinical.columns.difference(["CDR3", "Sample Type"])]
        localclinical = localclinical.loc[:,~localclinical.columns.duplicated()]
        localphysicochem = localphysicochem.loc[:,~localphysicochem.columns.duplicated()]
        
        
        localclinical_dfs = localclinical[localclinical.columns.difference(['OS_MONTHS', 'OS_STATUS'])]
        localclinical = localclinical[localclinical.columns.difference(['DFS_MONTHS', 'DFS_STATUS'])]
        
        localclinical[localclinical.columns.difference(['Filename', 'Sample', "Receptor"])] = localclinical[localclinical.columns.difference(['Filename', 'Sample', "Receptor"])].apply(pd.to_numeric, errors='coerce')
        localclinical = localclinical.dropna()
        localclinical = localclinical.set_index(["Receptor", "Filename"])
        
        localclinical_dfs[localclinical_dfs.columns.difference(['Filename', 'Sample', "Receptor"])] = localclinical_dfs[localclinical_dfs.columns.difference(['Filename', 'Sample', "Receptor"])].apply(pd.to_numeric, errors='coerce')
        localclinical_dfs = localclinical_dfs.dropna()
        localclinical_dfs = localclinical_dfs.set_index(["Receptor", "Filename"])
        
        tissue = localclinical.loc[localclinical["Sample"].str.contains("10")==False]
        tissue = tissue.loc[tissue["Sample"].str.contains("11")==False]
        tissue = tissue.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        tissue = tissue.reset_index()
        blood = localclinical.loc[localclinical["Sample"].str.contains("10")]
        blood = blood.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        blood = blood.reset_index()
        
        
        tissue_dfs = localclinical_dfs.loc[localclinical_dfs["Sample"].str.contains("10")==False]
        tissue_dfs = tissue_dfs.loc[tissue_dfs["Sample"].str.contains("11")==False]
        tissue_dfs = tissue_dfs.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        tissue_dfs = tissue_dfs.reset_index()
        blood_dfs = localclinical_dfs.loc[localclinical_dfs["Sample"].str.contains("10")]
        blood_dfs = blood_dfs.groupby(level = ["Receptor", "Filename"], axis=0).mean()
        blood_dfs = blood_dfs.reset_index()
        
        localphysicochem = localphysicochem.drop(["CDR3", "Sample", "Sample Type", "Receptor"], axis=1)
        self.km_physicochem_results = pd.DataFrame()
        
        
        for receptor in ["TR", "TRA|TRB", "TRG|TRD", "IGH|IGL", "IGH|IGK", "IG", "TRA", "TRB", "TRD", "TRG", "IGK", "IGL", "IGH"]:
            for ndiv in [2, 5]:
                for survivaltype in ["OS", "DFS"]:
                    for site in ["TISSUE", "BLOOD"]:
                        if survivaltype == "OS":
                            if site == "TISSUE":
                                try:
                                    results = kmcurve(localphysicochem, tissue, site, ndiv, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype, plot)
                                    self.km_physicochem_results = pd.concat([self.km_physicochem_results, results])
                                except:
                                    pass
                            if site == "BLOOD":
                                try:
                                    results = kmcurve(localphysicochem, blood, site, ndiv, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype, plot)
                                    self.km_physicochem_results = pd.concat([self.km_physicochem_results, results])
                                except:
                                    pass
                        if survivaltype == "DFS":
                            if site == "TISSUE":
                                try:
                                    results = kmcurve(localphysicochem, tissue_dfs, site, ndiv, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype, plot)
                                    self.km_physicochem_results = pd.concat([self.km_physicochem_results, results])
                                except:
                                    pass
                            if site == "BLOOD":
                                try:
                                    results = kmcurve(localphysicochem, blood_dfs, site, ndiv, receptor, "%s_MONTHS"%survivaltype, "%s_STATUS"%survivaltype, plot)
                                    self.km_physicochem_results = pd.concat([self.km_physicochem_results, results])
                                except:
                                    pass
        self.km_physicochem_results = self.km_physicochem_results.sort_values(by="logrank_pvalue")
        self.km_physicochem_results = self.km_physicochem_results.reset_index(drop=True)
        return self.km_physicochem_results
        
    def km_receptors(self, plot=False):
        def statscalc(df, df2, survivaltype, comparison_pvalues, plot):
            try:
                survivaltime = survivaltype+"_MONTHS"
                status = survivaltype+"_STATUS"
                df = df[[survivaltime, status]]
                df = df.apply(pd.to_numeric, errors='coerce')
                df = df.dropna()
                df = df.reset_index(drop=True)
                df2 = df2[[survivaltime, status]]
                df2 = df2.apply(pd.to_numeric, errors='coerce')
                df2 = df2.dropna()
                df2 = df2.reset_index(drop=True)
                logrank = logrank_test(df[survivaltime], df2[survivaltime], event_observed_A = df[status], event_observed_B = df2[status])
                n1 = len(df)
                n2 = len(df2)
                label = self.comparisonlist[len(comparison_pvalues)]
                label_A = label[:label.find("v")]
                label_B = label[label.find("v")+1:label.find("_")]
                median_res = self.medians(df[survivaltime], df2[survivaltime], event_observed_A = df[status], event_observed_B = df2[status], label_A = label_A, label_B = label_B)
                survival_avg1 = float(median_res[0])
                survival_avg2  = float(median_res[1])
                if plot == True:
                    fig = self.plotter(df[survivaltime], df2[survivaltime], event_observed_A = df[status], event_observed_B = df2[status], label_A = label_A, label_B = label_B)
                    fig.savefig(self.km_receptors_plots_directory+"/%s.pdf"%label)
                return logrank.p_value, logrank.test_statistic, survival_avg1, survival_avg2, n1, n2
            except:
                return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

                
        def quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n):
            print "Analyzing receptor KM curve: %s"%self.comparisonlist[len(comparison_pvalues)]
            comparison_pvalues.append(a[0])
            comparison_fvalues.append(a[1])
            group1_survival_avg.append(a[2])
            group2_survival_avg.append(a[3])
            group1_n.append(a[4])
            group2_n.append(a[5])
            return comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n
        
        localproductive = self.productive.set_index("Filename")
        localproductive = self.clinical[["Filename", "OS_STATUS", "OS_MONTHS", "DFS_MONTHS", "DFS_STATUS"]].join(localproductive, on='Filename')
        try:
            localproductive = localproductive.drop(["Reason"], axis = 1)
        except:
            pass
        localproductive = localproductive.dropna()
        localproductive = localproductive.reset_index()
        tissue = localproductive.loc[localproductive["Sample"].str.contains("10")==False]
        tissue = tissue.loc[tissue["Sample"].str.contains("11")==False]
        blood = localproductive.loc[localproductive["Sample"].str.contains("10")]
        self.comparisonlist = ["TRA|TRBvAR_os_tissue", 
        "TRA|TRBvAR_dfs_tissue", 
        "TRA|TRBvAR_os_blood", 
        "TRA|TRBvAR_dfs_blood", 
        "TRA|TRBvNoImmune_os_tissue", 
        "TRA|TRBvNoImmune_dfs_tissue", 
        "TRA|TRBvNoImmune_os_blood", 
        "TRA|TRBvNoImmune_dfs_blood", 
        "TRG|TRDvAR_os_tissue", 
        "TRG|TRDvAR_dfs_tissue", 
        "TRG|TRDvAR_os_blood", 
        "TRG|TRDvAR_dfs_blood", 
        "TRG|TRDvNoImmune_os_tissue", 
        "TRG|TRDvNoImmune_dfs_tissue", 
        "TRG|TRDvNoImmune_os_blood", 
        "TRG|TRDvNoImmune_dfs_blood", 
        "IGH|IGLvAR_os_tissue", 
        "IGH|IGLvAR_dfs_tissue", 
        "IGH|IGLvAR_os_blood", 
        "IGH|IGLvAR_dfs_blood", 
        "IGH|IGLvNoImmune_os_tissue", 
        "IGH|IGLvNoImmune_dfs_tissue", 
        "IGH|IGLvNoImmune_os_blood", 
        "IGH|IGLvNoImmune_dfs_blood", 
        "IGH|IGKvAR_os_tissue", 
        "IGH|IGKvAR_dfs_tissue", 
        "IGH|IGKvAR_os_blood", 
        "IGH|IGKvAR_dfs_blood", 
        "IGH|IGKvNoImmune_os_tissue", 
        "IGH|IGKvNoImmune_dfs_tissue", 
        "IGH|IGKvNoImmune_os_blood", 
        "IGH|IGKvNoImmune_dfs_blood", 
        "TR+IGvAR_os_tissue",
        "TRvAR_os_tissue",
        "TRAvAR_os_tissue",
        "TRBvAR_os_tissue",
        "TRDvAR_os_tissue",
        "TRGvAR_os_tissue",
        "IGvAR_os_tissue",
        "IGHvAR_os_tissue",
        "IGKvAR_os_tissue",
        "IGLvAR_os_tissue",
        "TR+IGvAR_dfs_tissue",
        "TRvAR_dfs_tissue",
        "TRAvAR_dfs_tissue",
        "TRBvAR_dfs_tissue",
        "TRDvAR_dfs_tissue",
        "TRGvAR_dfs_tissue",
        "IGvAR_dfs_tissue",
        "IGHvAR_dfs_tissue",
        "IGKvAR_dfs_tissue",
        "IGLvAR_dfs_tissue",
        "TRvNoImmune_os_tissue",
        "TRAvNoImmune_os_tissue",
        "TRBvNoImmune_os_tissue",
        "TRDvNoImmune_os_tissue",
        "TRGvNoImmune_os_tissue",
        "IGvNoImmune_os_tissue",
        "IGHvNoImmune_os_tissue",
        "IGKvNoImmune_os_tissue",
        "IGLvNoImmune_os_tissue",
        "TRvNoImmune_dfs_tissue",
        "TRAvNoImmune_dfs_tissue",
        "TRBvNoImmune_dfs_tissue",
        "TRDvNoImmune_dfs_tissue",
        "TRGvNoImmune_dfs_tissue",
        "IGvNoImmune_dfs_tissue",
        "IGHvNoImmune_dfs_tissue",
        "IGKvNoImmune_dfs_tissue",
        "IGLvNoImmune_dfs_tissue",
        "TRvIG_os_tissue",
        "TRvIG_dfs_tissue",
        "TR+IGvAR_os_blood",
        "TRvAR_os_blood",
        "TRAvAR_os_blood",
        "TRBvAR_os_blood",
        "TRDvAR_os_blood",
        "TRGvAR_os_blood",
        "IGvAR_os_blood",
        "IGHvAR_os_blood",
        "IGKvAR_os_blood",
        "IGLvAR_os_blood",
        "TR+IGvAR_dfs_blood",
        "TRvAR_dfs_blood",
        "TRAvAR_dfs_blood",
        "TRBvAR_dfs_blood",
        "TRDvAR_dfs_blood",
        "TRGvAR_dfs_blood",
        "IGvAR_dfs_blood",
        "IGHvAR_dfs_blood",
        "IGKvAR_dfs_blood",
        "IGLvAR_dfs_blood",
        "TRvNoImmune_os_blood",
        "TRAvNoImmune_os_blood",
        "TRBvNoImmune_os_blood",
        "TRDvNoImmune_os_blood",
        "TRGvNoImmune_os_blood",
        "IGvNoImmune_os_blood",
        "IGHvNoImmune_os_blood",
        "IGKvNoImmune_os_blood",
        "IGLvNoImmune_os_blood",
        "TRvNoImmune_dfs_blood",
        "TRAvNoImmune_dfs_blood",
        "TRBvNoImmune_dfs_blood",
        "TRDvNoImmune_dfs_blood",
        "TRGvNoImmune_dfs_blood",
        "IGvNoImmune_dfs_blood",
        "IGHvNoImmune_dfs_blood",
        "IGKvNoImmune_dfs_blood",
        "IGLvNoImmune_dfs_blood",
        "TRvIG_os_blood",
        "TRvIG_dfs_blood"]
        comparison_pvalues = []
        comparison_fvalues = []
        group1_survival_avg = []
        group2_survival_avg = []
        group1_n = []
        group2_n = []
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRA|TRB")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRA|TRB")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRA|TRB")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRA|TRB")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRA|TRB")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRG|TRD")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRG|TRD")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRG|TRD")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRG|TRD")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TRG|TRD")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGL")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGL")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGL")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGL")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGL")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGK")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGK")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGK")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGK")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IGH|IGK")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRA"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRB"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRD"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRG"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGH"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGK"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGL"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRA"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRB"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRD"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRG"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGH"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGK"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGL"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(tissue["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[self.clinical["Filename"].isin(tissue.loc[tissue["Receptor"].str.contains("IG")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRA"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRB"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRD"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRG"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGH"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGK"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGL"]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRA"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRB"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRD"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRG"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGH"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGK"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGL"]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRA"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRB"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRD"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="TRG"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGH"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)

        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGK"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"]=="IGL"]["Filename"])], self.clinical.loc[~self.clinical["Filename"].isin(blood["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], "OS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        a = statscalc(self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("TR")]["Filename"])], self.clinical.loc[self.clinical["Filename"].isin(blood.loc[blood["Receptor"].str.contains("IG")]["Filename"])], "DFS", comparison_pvalues, plot)
        comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n = quickappend(a, comparison_pvalues, comparison_fvalues, group1_survival_avg, group2_survival_avg, group1_n, group2_n)
        
        self.km_receptors_results = pd.DataFrame()
        self.km_receptors_results["comparison"] = self.comparisonlist
        self.km_receptors_results[["comparison", "survival_type", "site"]] = self.km_receptors_results['comparison'].str.split('_',expand=True)
        self.km_receptors_results["comparison_pvalues"] = comparison_pvalues
        self.km_receptors_results["comparison_fvalues"] = comparison_fvalues
        self.km_receptors_results["group1_survival_median"] = group1_survival_avg
        self.km_receptors_results["group2_survival_median"] = group2_survival_avg
        self.km_receptors_results["group1_group2_survival_difference"] = self.km_receptors_results["group1_survival_median"] - self.km_receptors_results["group2_survival_median"]
        self.km_receptors_results["group1_n"] = group1_n
        self.km_receptors_results["group2_n"] = group2_n
        self.km_receptors_results = self.km_receptors_results.sort_values(by="comparison_pvalues")
        self.km_receptors_results = self.km_receptors_results.reset_index(drop=True)
        return self.km_receptors_results
        