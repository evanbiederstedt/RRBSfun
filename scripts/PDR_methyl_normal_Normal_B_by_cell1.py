
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')



normal_B = glob.glob("RRBS_normal_B*")
newdf1 = pd.DataFrame()
for filename in normal_B:
    df = pd.read_table(filename)
    df['filename'] = str(filename)
    
    df = df.drop(['start', 'strand', 'avgWeightedEnt', 'CpGEntropy', 'tss', 'genes', 'exons', 'introns',
                  'promoter', 'cgi', 'geneDensity', 'ctcfUpstream',  'ctcfDownstream',
                  'ctcfDensity', 'geneDistalRegulatoryModules',  'vistaEnhancers', '3PrimeUTR', 'ctcfUpDistance', 'ctcfDownDistance',
                  '3PrimeUTRDistance', '5PrimeUTR', '5PrimeUTRDistance', 'firstExon',
                  'geneDistalRegulatoryModulesK562', 'geneDistalRegulatoryModulesK562Distance', 'hypoInHues64','hypoInHues64Distance',
                  'tssDistance', 'genesDistance', 'exonsDistance', 'intronsDistance', 'promoterDistance', 'cgiDistance',
                  'ctcf', 'ctcfDistance', 'geneDistalRegulatoryModulesDistance', 'vistaEnhancersDistance', 'firstExonDistance'], axis=1)

    chromosomes = ['chr2', 'chr5', 'chr11']
    df = df[(df["chr"].isin(chromosomes))]
                  
    num_bins = np.ceil(df["avgReadCpGs"].max()/1.25)
    df["avgReadCpGs_binned"] = pd.cut(df['avgReadCpGs'], num_bins, labels=False)
                  
    df["read_stack_ID"] = (df.avgReadCpGs_binned.shift(1) != df.avgReadCpGs_binned).astype(int).cumsum()
                  
    df["total_reads"] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1)
                  
    df = df.groupby(["read_stack_ID"])[["filename", "thisMeth", "thisUnmeth", "methReadCount", "unmethReadCount", "mixedReadCount", "total_reads"]].sum()
    df["filename"] = str(filename)
                  
    df = df[["filename", "thisMeth", "mixedReadCount", "total_reads"]].sum()
                  
    df["PDR_total"] = df["mixedReadCount"]/df["total_reads"]
                  
    df["methylation"] = df["thisMeth"]/(df["thisMeth"]+df["thisUnmeth"])          # corrected
                  
    newdf1 = newdf1.append(df, ignore_index=True)


newdf1 = newdf1[["filename", "methylation_total", "PDR_total", "thisMeth", "mixedReadCount", "total_reads"]]

# export as .csv
newdf1.to_csv("Meth_PDR_cell_RRBS_normal_B1.csv")