
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')


# 'tssDistance', 'genesDistance','exonsDistance', 'intronsDistance', 'promoterDistance', 'cgiDistance','ctcfDistance', 'ctcfUpDistance', 'ctcfDownDistance','geneDistalRegulatoryModulesDistance', 'vistaEnhancersDistance',  '3PrimeUTRDistance', '5PrimeUTRDistance',  'firstExonDistance','geneDistalRegulatoryModulesK562Distance', 'hypoInHues64Distance'


mcells = glob.glob("RRBS_NormalBCD19pCD27mcell*")
newdf1 = pd.DataFrame()
for filename in mcells:
    df = pd.read_table(filename)
    df['filename'] = str(filename)
    
    df = df.drop(['start', 'strand', 'avgWeightedEnt', 'CpGEntropy', 'tss', 'genes', 'exons', 'introns',
                  'promoter', 'cgi', 'geneDensity', 'ctcfUpstream',  'ctcfDownstream',
                  'ctcfDensity', 'geneDistalRegulatoryModules',  'vistaEnhancers', '3PrimeUTR', '5PrimeUTR', 'firstExon','geneDistalRegulatoryModulesK562', 'hypoInHues64','ctcf'], axis=1)
                  
    num_bins = np.ceil(df["avgReadCpGs"].max()/1.25)
    df["avgReadCpGs_binned"] = pd.cut(df['avgReadCpGs'], num_bins, labels=False)
                  
    df["read_stack_ID"] = (df.avgReadCpGs_binned.shift(1) != df.avgReadCpGs_binned).astype(int).cumsum()
                  
    df["total_reads"] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1)
                  
    df = df.groupby(["read_stack_ID"])[["filename", "thisMeth", "thisUnmeth", "methReadCount", "unmethReadCount", "mixedReadCount", "total_reads"]].sum()
    df["filename"] = str(filename)
                  
    df = df[["filename", "thisMeth", "thisUnmeth", "mixedReadCount", "total_reads"]].sum()
                  
    df["PDR_total"] = df["mixedReadCount"]/df["total_reads"]
                  
    df["methylation"] = df["thisMeth"]/(df["thisMeth"]+df["thisUnmeth"])          # corrected
                  
    newdf1 = newdf1.append(df, ignore_index=True)


newdf1 = newdf1[["filename", "methylation", "PDR_total", "thisMeth", "mixedReadCount", "total_reads", "ctcfUpDistance", "ctcfDownDistance", "3PrimeUTRDistance", "5PrimeUTRDistance", "geneDistalRegulatoryModulesK562Distance", "hypoInHues64Distance", "genesDistance", "exonsDistance", "intronsDistance", "promoterDistance", "cgiDistance", "ctcfDistance", "geneDistalRegulatoryModulesDistance", "vistaEnhancersDistance", "firstExonDistance", "tssDistance"]]

# export as .csv
newdf1.to_csv("Meth_PDR_cell_cd27mcell_genomic_regions.csv")