
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')

# set glob subdirectory via cell batch
cll_cells1 = glob.glob("RRBS_cw154*")

newdf1 = pd.DataFrame()
for filename in cll_cells1:
    df = pd.read_table(filename)
    df['filename'] = str(filename)
    df = df.drop(['start', 'strand', 'avgWeightedEnt', 'CpGEntropy', 'avgReadCpGs',
                  'tss', 'genes', 'exons', 'introns', 'promoter', 'cgi', 'geneDensity', 'ctcfUpstream',  'ctcfDownstream',
                  'ctcfDensity', 'geneDistalRegulatoryModules',  'vistaEnhancers', '3PrimeUTR', 'ctcfUpDistance', 'ctcfDownDistance',
                  '3PrimeUTRDistance', '5PrimeUTR', '5PrimeUTRDistance', 'firstExon',
                  'geneDistalRegulatoryModulesK562', 'geneDistalRegulatoryModulesK562Distance', 'hypoInHues64','hypoInHues64Distance',
                  'tssDistance', 'genesDistance', 'exonsDistance', 'intronsDistance', 'promoterDistance', 'cgiDistance',
                  'ctcf', 'ctcfDistance', 'geneDistalRegulatoryModulesDistance', 'vistaEnhancersDistance', 'firstExonDistance'], axis=1)
        
    chromosomes = ['chr2', 'chr5', 'chr11']
    df = df[(df['chr'].isin(chromosomes))]
                  
    df["total_reads"] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1)
                  
    df = df.sum()
    # Methylation == # methylated reads per CpG site / total # of reads per CpG site
    df["percent_methylation"] = df["methReadCount"]/df["total_reads"]
    newdf1 = newdf1.append(df, ignore_index=True)

newdf1 = newdf1.drop(['chr'], axis=1)

newdf1 = newdf1[['filename', 'percent_methylation', 'thisMeth', 'thisUnmeth', 'methReadCount', 'unmethReadCount', 'mixedReadCount', 'total_reads']]

# export as .csv
newdf1.to_csv('CLL_cells_methylation_chr2_5_11A.csv')
