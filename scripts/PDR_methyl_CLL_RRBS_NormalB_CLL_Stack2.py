
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')


# set glob subdirectory via cell batch
cll_cells2 = glob.glob("RRBS_NormalB_CLL*")


for filename in cll_cells2:
    df = pd.read_table(filename)
    
    df = df.drop(['chr', 'start', 'strand', 'avgWeightedEnt', 'CpGEntropy', 'tss', 'tssDistance', 'genes', 'exons',
                  'introns', 'promoter', 'cgi',
                  'geneDensity', 'ctcfUpstream',  'ctcfDownstream','ctcfDensity', 'geneDistalRegulatoryModules',
                  'vistaEnhancers', '3PrimeUTR', 'ctcfUpDistance', 'ctcfDownDistance','3PrimeUTRDistance',
                  '5PrimeUTR', '5PrimeUTRDistance', 'firstExon','geneDistalRegulatoryModulesK562',
                  'geneDistalRegulatoryModulesK562Distance', 'hypoInHues64','hypoInHues64Distance',
                  'genesDistance', 'exonsDistance', 'intronsDistance', 'promoterDistance', 'cgiDistance',
                  'ctcf', 'ctcfDistance', 'geneDistalRegulatoryModulesDistance', 'vistaEnhancersDistance', 'firstExonDistance'], axis=1)



    num_bins2 = np.ceil(df['avgReadCpGs'].max()/1.25)
    df['avgReadCpGs_binned'] = pd.cut(df['avgReadCpGs'], num_bins2, labels=False)

    df['read_stack_ID'] = (df.avgReadCpGs_binned.shift(1) != df.avgReadCpGs_binned).astype(int).cumsum()

    df['total_reads'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1)


    df['avgReadCpGs'] = df['avgReadCpGs'].values.round(decimals=0)

    df1 = df.groupby(['read_stack_ID', 'avgReadCpGs'])[['thisMeth', 'thisUnmeth', 'methReadCount', 'unmethReadCount', 'mixedReadCount', 'total_reads']].sum()

    df1.reset_index(inplace=True)


    df1["percent_methylationPerStack"] = df1["thisMeth"]/df1["total_reads"]
    df1["PDR_per_stack"] = df1["mixedReadCount"]/df1["total_reads"]
    
    df1.to_csv(str("stacked_") + str(filename) +str(".csv"))
                  
