import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')

pcell = glob.glob("RRBS_NormalBCD19pCD27pcell*")

newdf1 = pd.DataFrame()
for filename in pcell:
    df = pd.read_table(filename)
    df['filename'] = str(filename)
    
    df = df.drop(['chr', 'start', 'strand', 'avgWeightedEnt','CpGEntropy', 'tss', 'genes', 'exons', 'introns', 'promoter',
       'cgi', 'geneDensity', 'ctcf','ctcfUpstream', 'ctcfDownstream', 'ctcfDensity', 'geneDistalRegulatoryModules',
       'vistaEnhancers','3PrimeUTR',  '5PrimeUTR','firstExon', 'geneDistalRegulatoryModulesK562','hypoInHues64'], axis=1)
    
    num_bins2 = np.ceil(df['avgReadCpGs'].max()/1.25)
    
    df['avgReadCpGs_binned'] = pd.cut(df['avgReadCpGs'], num_bins2, labels=False)
    
    df['read_stack_ID'] = (df.avgReadCpGs_binned.shift(1) != df.avgReadCpGs_binned).astype(int).cumsum()
    df['total_reads'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1)
    df['avgReadCpGs'] = df['avgReadCpGs'].values.round(decimals=0)
    
    intergenic_set = ( (df['genesDistance']>0) & (df['exonsDistance']>0) & (df['intronsDistance']>0) & (df['promoterDistance']>0) & (df['cgiDistance']>0) & (df['ctcfDistance']>0) & 
 (df['ctcfUpDistance']>0) & (df['ctcfDownDistance']>0) & (df['geneDistalRegulatoryModulesDistance']>0) & (df['vistaEnhancersDistance']>0) & 
 (df['3PrimeUTRDistance']>0) & (df['5PrimeUTRDistance']>0) & (df['firstExonDistance']>0) & (df['geneDistalRegulatoryModulesK562Distance']>0) & (df['hypoInHues64Distance']>0) )
    
    df['totreads_tssDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['tssDistance']<0, 0)
    df['totreads_genesDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['genesDistance']<0, 0)
    df['totreads_exonsDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['exonsDistance']<0, 0)
    df['totreads_intronsDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['intronsDistance']<0, 0)
    df['totreads_promoterDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['promoterDistance']<0, 0)
    df['totreads_cgiDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['cgiDistance']<0, 0)
    df['totreads_ctcfDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['ctcfDistance']<0, 0)
    df['totreads_ctcfUpDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['ctcfUpDistance']<0, 0)
    df['totreads_ctcfDownDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['ctcfDownDistance']<0, 0)
    df['totreads_geneDistalRegulatoryModulesDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['geneDistalRegulatoryModulesDistance']<0, 0)
    df['totreads_vistaEnhancersDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['vistaEnhancersDistance']<0, 0)
    df['totreads_3PrimeUTRDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['3PrimeUTRDistance']<0, 0)
    df['totreads_5PrimeUTRDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['5PrimeUTRDistance']<0, 0)
    df['totreads_firstExonDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['firstExonDistance']<0, 0)
    df['totreads_geneDistalRegulatoryModulesK562Distance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['geneDistalRegulatoryModulesK562Distance']<0, 0)
    df['totreads_hypoInHues64Distance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['hypoInHues64Distance']<0, 0)
    df['totreads_intergenic'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(intergenic_set, 0)
    #
    # shore (+1 to +2000....on either side) 
    # shelf (+2000 to +4000....on either side)
    #
    df['totreads_shore'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where((df['cgiDistance'].between(0, 2000, inclusive=False)), 0)  # 0 < cgidistance < 2000
    df['totreads_shelf'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(((df['cgiDistance'] >= 2000) & (df['cgiDistance'] < 4000)), 0) # 2000 =< cgidistance < 4000

    
    df['mixedReads_tssDistance'] = np.where(df['tssDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_genesDistance'] = np.where(df['genesDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_exonsDistance'] = np.where(df['exonsDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_intronsDistance'] = np.where(df['intronsDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_promoterDistance'] = np.where(df['promoterDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_cgiDistance'] = np.where(df['cgiDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_ctcfDistance'] = np.where(df['ctcfDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_ctcfUpDistance'] = np.where(df['ctcfUpDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_ctcfDownDistance'] = np.where(df['ctcfDownDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_geneDistalRegulatoryModulesDistance'] = np.where(df['geneDistalRegulatoryModulesDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_vistaEnhancersDistance'] = np.where(df['vistaEnhancersDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_3PrimeUTRDistance'] = np.where(df['3PrimeUTRDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_5PrimeUTRDistance'] = np.where(df['5PrimeUTRDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_firstExonDistance'] = np.where(df['firstExonDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_geneDistalRegulatoryModulesK562Distance'] = np.where(df['geneDistalRegulatoryModulesK562Distance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_hypoInHues64Distance'] = np.where(df['hypoInHues64Distance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_intergenic'] = np.where(intergenic_set, df['mixedReadCount'], 0)
    df['mixedReads_shore'] = np.where((df['cgiDistance'].between(0, 2000, inclusive=False)), df['mixedReadCount'], 0)
    df['mixedReads_shelf'] = np.where(((df['cgiDistance'] >= 2000) & (df['cgiDistance'] < 4000)), df['mixedReadCount'], 0)

    
    df['thisMeth_tssDistance'] = np.where(df['tssDistance']<0, df['thisMeth'], 0)    
    df['thisMeth_genesDistance'] = np.where(df['genesDistance']<0, df['thisMeth'], 0)
    df['thisMeth_exonsDistance'] = np.where(df['exonsDistance']<0, df['thisMeth'], 0)
    df['thisMeth_intronsDistance'] = np.where(df['intronsDistance']<0, df['thisMeth'], 0)
    df['thisMeth_promoterDistance'] = np.where(df['promoterDistance']<0, df['thisMeth'], 0)
    df['thisMeth_cgiDistance'] = np.where(df['cgiDistance']<0, df['thisMeth'], 0)
    df['thisMeth_ctcfDistance'] = np.where(df['ctcfDistance']<0, df['thisMeth'], 0)
    df['thisMeth_ctcfUpDistance'] = np.where(df['ctcfUpDistance']<0, df['thisMeth'], 0)
    df['thisMeth_ctcfDownDistance'] = np.where(df['ctcfDownDistance']<0, df['thisMeth'], 0)
    df['thisMeth_geneDistalRegulatoryModulesDistance'] = np.where(df['geneDistalRegulatoryModulesDistance'] <0, df['thisMeth'], 0)
    df['thisMeth_vistaEnhancersDistance'] = np.where(df['vistaEnhancersDistance'] <0, df['thisMeth'], 0)
    df['thisMeth_3PrimeUTRDistance'] = np.where(df['3PrimeUTRDistance'] <0, df['thisMeth'], 0)
    df['thisMeth_5PrimeUTRDistance'] = np.where(df['5PrimeUTRDistance'] <0, df['thisMeth'], 0)
    df['thisMeth_firstExonDistance'] = np.where(df['firstExonDistance'] <0, df['thisMeth'], 0)
    df['thisMeth_geneDistalRegulatoryModulesK562Distance'] = np.where(df['geneDistalRegulatoryModulesK562Distance'] <0, df['thisMeth'], 0)
    df['thisMeth_hypoInHues64Distance'] = np.where(df['hypoInHues64Distance'] <0, df['thisMeth'], 0)
    df['thisMeth_intergenic'] = np.where(intergenic_set, df['thisMeth'], 0)
    df['thisMeth_shore'] = np.where((df['cgiDistance'].between(0, 2000, inclusive=False)), df['thisMeth'], 0)
    df['thisMeth_shelf'] = np.where(((df['cgiDistance'] >= 2000) & (df['cgiDistance'] < 4000)), df['thisMeth'], 0)
    
    
    df['thisUnmeth_tssDistance'] = np.where(df['tssDistance']<0, df['thisUnmeth'], 0)   
    df['thisUnmeth_genesDistance'] = np.where(df['genesDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_exonsDistance'] = np.where(df['exonsDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_intronsDistance'] = np.where(df['intronsDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_promoterDistance'] = np.where(df['promoterDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_cgiDistance'] = np.where(df['cgiDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_ctcfDistance'] = np.where(df['ctcfDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_ctcfUpDistance'] = np.where(df['ctcfUpDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_ctcfDownDistance'] = np.where(df['ctcfDownDistance']<0, df['thisUnmeth'], 0)
    df['thisUnmeth_geneDistalRegulatoryModulesDistance'] = np.where(df['geneDistalRegulatoryModulesDistance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_vistaEnhancersDistance'] = np.where(df['vistaEnhancersDistance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_3PrimeUTRDistance'] = np.where(df['3PrimeUTRDistance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_5PrimeUTRDistance'] = np.where(df['5PrimeUTRDistance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_firstExonDistance'] = np.where(df['firstExonDistance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_geneDistalRegulatoryModulesK562Distance'] = np.where(df['geneDistalRegulatoryModulesK562Distance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_hypoInHues64Distance'] = np.where(df['hypoInHues64Distance'] <0, df['thisUnmeth'], 0)
    df['thisUnmeth_intergenic'] = np.where(intergenic_set, df['thisUnmeth'], 0)
    df['thisUnmeth_shore'] = np.where((df['cgiDistance'].between(0, 2000, inclusive=False)), df['thisUnmeth'], 0)
    df['thisUnmeth_shelf'] = np.where(((df['cgiDistance'] >= 2000) & (df['cgiDistance'] < 4000)), df['thisUnmeth'], 0)
    
    df = df.groupby(['read_stack_ID'])[["filename", 'thisMeth', 'thisUnmeth', 'methReadCount', 'unmethReadCount', 'mixedReadCount', 
        'total_reads','totreads_tssDistance', 'totreads_genesDistance', 'totreads_exonsDistance', 'totreads_intronsDistance','totreads_promoterDistance', 'totreads_cgiDistance',
       'totreads_ctcfDistance', 'totreads_ctcfUpDistance','totreads_ctcfDownDistance','totreads_geneDistalRegulatoryModulesDistance',
       'totreads_vistaEnhancersDistance', 'totreads_3PrimeUTRDistance','totreads_5PrimeUTRDistance', 'totreads_firstExonDistance',
       'totreads_geneDistalRegulatoryModulesK562Distance','totreads_hypoInHues64Distance', 'totreads_intergenic','totreads_shore', 
       'totreads_shelf', 'mixedReads_tssDistance', 'mixedReads_genesDistance','mixedReads_exonsDistance', 'mixedReads_intronsDistance','mixedReads_promoterDistance', 
       'mixedReads_cgiDistance','mixedReads_ctcfDistance', 'mixedReads_ctcfUpDistance','mixedReads_ctcfDownDistance',
       'mixedReads_geneDistalRegulatoryModulesDistance','mixedReads_vistaEnhancersDistance', 'mixedReads_3PrimeUTRDistance',
       'mixedReads_5PrimeUTRDistance', 'mixedReads_firstExonDistance','mixedReads_geneDistalRegulatoryModulesK562Distance','mixedReads_hypoInHues64Distance', 'mixedReads_intergenic',
       'mixedReads_shore', 'mixedReads_shelf', 'thisMeth_tssDistance','thisMeth_genesDistance','thisMeth_exonsDistance', 'thisMeth_intronsDistance',
       'thisMeth_promoterDistance', 'thisMeth_cgiDistance','thisMeth_ctcfDistance', 'thisMeth_ctcfUpDistance','thisMeth_ctcfDownDistance',
       'thisMeth_geneDistalRegulatoryModulesDistance','thisMeth_vistaEnhancersDistance', 'thisMeth_3PrimeUTRDistance',
       'thisMeth_5PrimeUTRDistance', 'thisMeth_firstExonDistance','thisMeth_geneDistalRegulatoryModulesK562Distance',
       'thisMeth_hypoInHues64Distance', 'thisMeth_intergenic','thisMeth_shore', 'thisMeth_shelf', 'thisUnmeth_tssDistance','thisUnmeth_genesDistance',
       'thisUnmeth_exonsDistance', 'thisUnmeth_intronsDistance','thisUnmeth_promoterDistance', 'thisUnmeth_cgiDistance',
       'thisUnmeth_ctcfDistance', 'thisUnmeth_ctcfUpDistance','thisUnmeth_ctcfDownDistance','thisUnmeth_geneDistalRegulatoryModulesDistance',
       'thisUnmeth_vistaEnhancersDistance', 'thisUnmeth_3PrimeUTRDistance','thisUnmeth_5PrimeUTRDistance', 'thisUnmeth_firstExonDistance',
       'thisUnmeth_geneDistalRegulatoryModulesK562Distance','thisUnmeth_hypoInHues64Distance', 'thisUnmeth_intergenic','thisUnmeth_shore', 'thisUnmeth_shelf']].sum()
    
    df = df[['thisMeth', 'thisUnmeth', 'methReadCount', 'unmethReadCount', 'mixedReadCount', 
        'total_reads', 'totreads_tssDistance', 'totreads_genesDistance', 'totreads_exonsDistance', 'totreads_intronsDistance','totreads_promoterDistance', 'totreads_cgiDistance',
       'totreads_ctcfDistance', 'totreads_ctcfUpDistance','totreads_ctcfDownDistance','totreads_geneDistalRegulatoryModulesDistance',
       'totreads_vistaEnhancersDistance', 'totreads_3PrimeUTRDistance','totreads_5PrimeUTRDistance', 'totreads_firstExonDistance',
       'totreads_geneDistalRegulatoryModulesK562Distance','totreads_hypoInHues64Distance', 'totreads_intergenic','totreads_shore', 
       'totreads_shelf', 'mixedReads_tssDistance', 'mixedReads_genesDistance','mixedReads_exonsDistance', 'mixedReads_intronsDistance','mixedReads_promoterDistance', 
       'mixedReads_cgiDistance','mixedReads_ctcfDistance', 'mixedReads_ctcfUpDistance','mixedReads_ctcfDownDistance',
       'mixedReads_geneDistalRegulatoryModulesDistance','mixedReads_vistaEnhancersDistance', 'mixedReads_3PrimeUTRDistance',
       'mixedReads_5PrimeUTRDistance', 'mixedReads_firstExonDistance','mixedReads_geneDistalRegulatoryModulesK562Distance','mixedReads_hypoInHues64Distance', 'mixedReads_intergenic',
       'mixedReads_shore', 'mixedReads_shelf', 'thisMeth_tssDistance', 'thisMeth_genesDistance','thisMeth_exonsDistance', 'thisMeth_intronsDistance',
       'thisMeth_promoterDistance', 'thisMeth_cgiDistance','thisMeth_ctcfDistance', 'thisMeth_ctcfUpDistance','thisMeth_ctcfDownDistance',
       'thisMeth_geneDistalRegulatoryModulesDistance','thisMeth_vistaEnhancersDistance', 'thisMeth_3PrimeUTRDistance',
       'thisMeth_5PrimeUTRDistance', 'thisMeth_firstExonDistance','thisMeth_geneDistalRegulatoryModulesK562Distance',
       'thisMeth_hypoInHues64Distance', 'thisMeth_intergenic','thisMeth_shore', 'thisMeth_shelf', 'thisUnmeth_tssDistance','thisUnmeth_genesDistance',
       'thisUnmeth_exonsDistance', 'thisUnmeth_intronsDistance','thisUnmeth_promoterDistance', 'thisUnmeth_cgiDistance',
       'thisUnmeth_ctcfDistance', 'thisUnmeth_ctcfUpDistance','thisUnmeth_ctcfDownDistance','thisUnmeth_geneDistalRegulatoryModulesDistance',
       'thisUnmeth_vistaEnhancersDistance', 'thisUnmeth_3PrimeUTRDistance','thisUnmeth_5PrimeUTRDistance', 'thisUnmeth_firstExonDistance',
       'thisUnmeth_geneDistalRegulatoryModulesK562Distance','thisUnmeth_hypoInHues64Distance', 'thisUnmeth_intergenic','thisUnmeth_shore', 'thisUnmeth_shelf']].sum()
    
    df["filename"] = str(filename)


    try:
        df["methylation_tssDistance"] = df["thisMeth_tssDistance"]/(df["thisMeth_tssDistance"]+df["thisUnmeth_tssDistance"])  
    except ZeroDivisionError:
        df["methylation_tssDistance"] = 0                                            
    try:
        df["methylation_genesDistance"] = df["thisMeth_genesDistance"]/(df["thisMeth_genesDistance"]+df["thisUnmeth_genesDistance"])  
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0
    try:
        df["methylation_exonsDistance"] = df["thisMeth_exonsDistance"]/(df["thisMeth_exonsDistance"]+df["thisUnmeth_exonsDistance"])  
    except ZeroDivisionError:
        df["methylation_exonsDistance"] = 0    
    try:   
        df["methylation_intronsDistance"] = df["thisMeth_intronsDistance"]/(df["thisMeth_intronsDistance"]+df["thisUnmeth_intronsDistance"])  
    except ZeroDivisionError:
        df["methylation_intronsDistance"] = 0            
    try:   
        df["methylation_promoterDistance"] = df["thisMeth_promoterDistance"]/(df["thisMeth_promoterDistance"]+df["thisUnmeth_promoterDistance"])  
    except ZeroDivisionError:
        df["methylation_promoterDistance"] = 0            
    try:       
        df["methylation_cgiDistance"] = df["thisMeth_cgiDistance"]/(df["thisMeth_cgiDistance"]+df["thisUnmeth_cgiDistance"])  
    except ZeroDivisionError:
        df["methylation_cgiDistance"] = 0            
    try:      
        df["methylation_ctcfDistance"] = df["thisMeth_ctcfDistance"]/(df["thisMeth_ctcfDistance"]+df["thisUnmeth_ctcfDistance"])  
    except ZeroDivisionError:
        df["methylation_ctcfDistance"] = 0         
    try:       
        df["methylation_ctcfUpDistance"] = df["thisMeth_ctcfUpDistance"]/(df["thisMeth_ctcfUpDistance"]+df["thisUnmeth_ctcfUpDistance"])  
    except ZeroDivisionError:
        df["methylation_ctcfUpDistance"] = 0         
    try:        
        df["methylation_ctcfDownDistance"] = df["thisMeth_ctcfDownDistance"]/(df["thisMeth_ctcfDownDistance"]+df["thisUnmeth_ctcfDownDistance"])  
    except ZeroDivisionError:
        df["methylation_ctcfDownDistance"] = 0         
    try:        
        df["methylation_geneDistalRegulatoryModulesDistance"] = df["thisMeth_geneDistalRegulatoryModulesDistance"]/(df["thisMeth_geneDistalRegulatoryModulesDistance"]+df["thisUnmeth_geneDistalRegulatoryModulesDistance"])  
    except ZeroDivisionError:
        df["methylation_geneDistalRegulatoryModulesDistance"] = 0         
    try:        
        df["methylation_vistaEnhancersDistance"] = df["thisMeth_vistaEnhancersDistance"]/(df["thisMeth_vistaEnhancersDistance"]+df["thisUnmeth_vistaEnhancersDistance"])  
    except ZeroDivisionError:
        df["methylation_vistaEnhancersDistance"] = 0         
    try:        
        df["methylation_3PrimeUTRDistance"] = df["thisMeth_3PrimeUTRDistance"]/(df["thisMeth_3PrimeUTRDistance"]+df["thisUnmeth_3PrimeUTRDistance"])  
    except ZeroDivisionError:
        df["methylation_3PrimeUTRDistance"] = 0         
    try:        
        df["methylation_5PrimeUTRDistance"] = df["thisMeth_5PrimeUTRDistance"]/(df["thisMeth_5PrimeUTRDistance"]+df["thisUnmeth_5PrimeUTRDistance"])  
    except ZeroDivisionError:
        df["methylation_PrimeUTRDistance"] = 0         
    try:        
        df["methylation_firstExonDistance"] = df["thisMeth_firstExonDistance"]/(df["thisMeth_firstExonDistance"]+df["thisUnmeth_firstExonDistance"])  
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0         
    try:        
        df["methylation_geneDistalRegulatoryModulesK562Distance"] = df["thisMeth_geneDistalRegulatoryModulesK562Distance"]/(df["thisMeth_geneDistalRegulatoryModulesK562Distance"]+df["thisUnmeth_geneDistalRegulatoryModulesK562Distance"])  
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0             
    try:        
        df["methylation_hypoInHues64Distance"] = df["thisMeth_hypoInHues64Distance"]/(df["thisMeth_hypoInHues64Distance"]+df["thisUnmeth_hypoInHues64Distance"])  
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0         
    try:        
        df["methylation_intergenic"] = df["thisMeth_intergenic"]/(df["thisMeth_intergenic"]+df["thisUnmeth_intergenic"]) 
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0         
    try:        
        df["methylation_shore"] = df["thisMeth_shore"]/(df["thisMeth_shore"]+df["thisUnmeth_shore"])  
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0         
    try:        
        df["methylation_shelf"] = df["thisMeth_shelf"]/(df["thisMeth_shelf"]+df["thisUnmeth_shelf"])  
    except ZeroDivisionError:
        df["methylation_genesDistance"] = 0 
        
        
        
        
        
        
        
    try:    
        df["PDR_tssDistance"] = df["mixedReads_genesDistance"]/df["totreads_tssDistance"]  
    except ZeroDivisionError:
        df["PDR_tssDistance"] = 0  
    try: 
        df["PDR_genesDistance"] = df["mixedReads_genesDistance"]/df["totreads_genesDistance"]
    except ZeroDivisionError:
        df["PDR_genesDistance"] = 0         
    try: 
        df["PDR_exonsDistance"] = df["mixedReads_exonsDistance"]/df["totreads_exonsDistance"]
    except ZeroDivisionError:
        df["PDR_exonsDistance"] = 0         
    try:
        df["PDR_intronsDistance"] = df["mixedReads_intronsDistance"]/df["totreads_intronsDistance"]
    except ZeroDivisionError:
        df["PDR_intronsDistance"] = 0 
    try:    
        df["PDR_promoterDistance"] = df["mixedReads_promoterDistance"]/df["totreads_promoterDistance"]
    except ZeroDivisionError:
        df["PDR_promoterDistance"] = 0          
    try:
        df["PDR_cgiDistance"] = df["mixedReads_cgiDistance"]/df["totreads_cgiDistance"]
    except ZeroDivisionError:
        df["PDR_cgiDistance"] = 0          
    try: 
        df["PDR_ctcfDistance"] = df["mixedReads_ctcfDistance"]/df["totreads_ctcfDistance"]
    except ZeroDivisionError:
        df["PDR_ctcfDistance"] = 0          
    try: 
        df["PDR_ctcfUpDistance"] = df["mixedReads_ctcfUpDistance"]/df["totreads_ctcfUpDistance"]
    except ZeroDivisionError:
        df["PDR_ctcfUpDistance"] = 0          
    try: 
        df["PDR_ctcfDownDistance"] = df["mixedReads_ctcfDownDistance"]/df["totreads_ctcfDownDistance"]
    except ZeroDivisionError:
        df["PDR_ctcfDownDistance"] = 0          
    try: 
        df["PDR_geneDistalRegulatoryModulesDistance"] = df["mixedReads_geneDistalRegulatoryModulesDistance"]/df["totreads_geneDistalRegulatoryModulesDistance"]
    except ZeroDivisionError:
        df["PDR_geneDistalRegulatoryModulesDistance"] = 0          
    try:
        df["PDR_vistaEnhancersDistance"] = df["mixedReads_vistaEnhancersDistance"]/df["totreads_vistaEnhancersDistance"]
    except ZeroDivisionError:
        df["PDR_vistaEnhancersDistance"] = 0          
    try:
        df["PDR_3PrimeUTRDistance"] = df["mixedReads_3PrimeUTRDistance"]/df["totreads_3PrimeUTRDistance"]
    except ZeroDivisionError:
        df["PDR_3PrimeUTRDistance"] = 0          
    try:
        df["PDR_5PrimeUTRDistance"] = df["mixedReads_5PrimeUTRDistance"]/df["totreads_5PrimeUTRDistance"]
    except ZeroDivisionError:
        df["PDR_5PrimeUTRDistance"] = 0          
    try: 
        df["PDR_firstExonDistance"] = df["mixedReads_firstExonDistance"]/df["totreads_firstExonDistance"]
    except ZeroDivisionError:
        df["PDR_firstExonDistance"] = 0          
    try: 
        df["PDR_geneDistalRegulatoryModulesK562Distance"] = df["mixedReads_geneDistalRegulatoryModulesK562Distance"]/df["totreads_geneDistalRegulatoryModulesK562Distance"]
    except ZeroDivisionError:
        df["PDR_geneDistalRegulatoryModulesK562Distance"] = 0          
    try: 
        df["PDR_hypoInHues64Distance"] = df["mixedReads_hypoInHues64Distance"]/df["totreads_hypoInHues64Distance"]
    except ZeroDivisionError:
        df["PDR_hypoInHues64Distance"] = 0          
    try:
        df["PDR_intergenic"] = df["mixedReads_intergenic"]/df["totreads_intergenic"]
    except ZeroDivisionError:
        df["PDR_intergenic"] = 0          
    try:
        df["PDR_shore"] = df["mixedReads_shore"]/df["totreads_shore"]
    except ZeroDivisionError:
        df["PDR_shore"] = 0          
    try: 
        df["PDR_shelf"] = df["mixedReads_shelf"]/df["totreads_shelf"]
    except ZeroDivisionError:
        df["PDR_shelf"] = 0  
        
    df["filename"] = str(filename)
        
    newdf1 = newdf1.append(df, ignore_index=True)



newdf1 = newdf1[["filename", "methylation_tssDistance","methylation_genesDistance","methylation_exonsDistance","methylation_intronsDistance",
                 "methylation_promoterDistance","methylation_cgiDistance","methylation_ctcfDistance","methylation_ctcfUpDistance","methylation_ctcfDownDistance",
             "methylation_geneDistalRegulatoryModulesDistance","methylation_vistaEnhancersDistance","methylation_3PrimeUTRDistance","methylation_5PrimeUTRDistance", "methylation_firstExonDistance","methylation_geneDistalRegulatoryModulesK562Distance","methylation_hypoInHues64Distance","methylation_intergenic",
                 "methylation_shore","methylation_shelf","PDR_tssDistance","PDR_genesDistance","PDR_exonsDistance","PDR_intronsDistance", "PDR_promoterDistance",
                 "PDR_cgiDistance", "PDR_ctcfDistance","PDR_ctcfUpDistance","PDR_ctcfDownDistance","PDR_geneDistalRegulatoryModulesDistance",
                 "PDR_vistaEnhancersDistance","PDR_3PrimeUTRDistance", "PDR_5PrimeUTRDistance","PDR_firstExonDistance","PDR_geneDistalRegulatoryModulesK562Distance",
                 "PDR_hypoInHues64Distance","PDR_intergenic","PDR_shore","PDR_shelf"]]

newdf1.to_csv("pcell_regions.csv", index=False)



