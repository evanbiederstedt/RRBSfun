
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')


import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')


# set glob subdirectory via cell batch
cll_cells3 = glob.glob("RRBS_trito_pool*")


newdf1 = pd.DataFrame()
for filename in cll_cells3:
    df = pd.read_table(filename)
    df['filename'] = str(filename)
    
    df = df.drop(['chr', 'start', 'strand', 'avgWeightedEnt', 'CpGEntropy', 'tss', 'tssDistance', 'genes', 'exons',
                  'introns', 'promoter', 'cgi','geneDensity', 'ctcfUpstream',  'ctcfDownstream','ctcfDensity',
                  'geneDistalRegulatoryModules','vistaEnhancers', '3PrimeUTR', 'ctcfUpDistance', 'ctcfDownDistance','3PrimeUTRDistance',
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
                  
    df1["methylation_per_stack"] = df1["thisMeth"]/(df1["thisMeth"]+df1["thisUnmeth"])          # corrected
    df1["PDR_per_stack_unweighted"] = (df1["mixedReadCount"]/df1["total_reads"])
                  
    df1["filename"] = str(filename)
                  
    df2 = df1[["filename", "methylation_per_stack", "PDR_per_stack_unweighted"]].sum()
                  
    df2["methylation_unweighted"] = df2["methylation_per_stack"]/(df1.read_stack_ID.max())
    df2["PDR_unweighted"] = df2["PDR_per_stack_unweighted"]/(df1.read_stack_ID.max())
                  
                  
    newdf1 = newdf1.append(df2, ignore_index=True)



newdf1 = newdf1[["filename", "PDR_unweighted", "methylation_unweighted"]]

# export as .csv
newdf1.to_csv("unweighted_PDR_methyl_RRBS_trito_pool.csv")


