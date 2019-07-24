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


df = pd.read_excel("")

#df = df.loc[df["Barcode"] != "TCGA-R8-A6YH"]


count = 0


countlist = []
p_value_list = []
while count < 1000:
    count += 1
    df = df.dropna()
    
    df1 = df.sample(129)
    df2 = df.loc[~df["Barcode"].isin(df1["Barcode"])]

    
    logrank = logrank_test(df1["OS_MONTHS"], df2["OS_MONTHS"], event_observed_A = df1["OS_STATUS"], event_observed_B = df2["OS_STATUS"])
    
    print count
    countlist.append(count)
    p_value_list.append(logrank.p_value)
    
simulation = pd.DataFrame()
simulation["Simulation Number"] = countlist
simulation["P Value"] = p_value_list
simulation.to_csv("")
