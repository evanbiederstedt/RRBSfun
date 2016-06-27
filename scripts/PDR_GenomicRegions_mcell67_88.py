
import glob
import pandas as pd
import numpy as np

import os
os.chdir('/Users/evanbiederstedt/Downloads/annoOld_files')

# set glob subdirectory via cell batch
normalB_cellbatch_mcell = glob.glob("RRBS_NormalBCD19pCD27mcell67_88*.anno")


newdf1 = pd.DataFrame()
for filename in normalB_cellbatch_mcell:
    df = pd.read_table(filename)
    df['filename'] = str(filename)
    df = df.drop(['chr', 'start', 'strand', 'thisMeth', 'thisUnmeth', 'avgWeightedEnt', 'CpGEntropy',
                  'avgReadCpGs', 'tss', 'tssDistance', 'genes', 'exons', 'introns', 'promoter', 'cgi',
                  'geneDensity', 'ctcfUpstream',  'ctcfDownstream','ctcfDensity', 'geneDistalRegulatoryModules',
                  'vistaEnhancers', '3PrimeUTR', 'ctcfUpDistance', 'ctcfDownDistance','3PrimeUTRDistance',
                  '5PrimeUTR', '5PrimeUTRDistance', 'firstExon','geneDistalRegulatoryModulesK562',
                  'geneDistalRegulatoryModulesK562Distance', 'hypoInHues64','hypoInHues64Distance'], axis=1)

    df['total_reads'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1)
    
    df['totreads_genesDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['genesDistance']<0, 0)
    df['totreads_exonsDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['exonsDistance']<0, 0)
    df['totreads_intronsDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['intronsDistance']<0, 0)
    df['totreads_promoterDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['promoterDistance']<0, 0)
    df['totreads_cgiDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['cgiDistance']<0, 0)
    df['totreads_ctcfDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['ctcfDistance']<0, 0)
    df['totreads_geneDistalRegulatoryModulesDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['geneDistalRegulatoryModulesDistance']<0, 0)
    df['totreads_firstExonDistance'] = df[["methReadCount", "unmethReadCount", "mixedReadCount"]].sum(axis=1).where(df['firstExonDistance']<0, 0)
                                    
                                    
    df['mixedReads_genesDistance'] = np.where(df['genesDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_exonsDistance'] = np.where(df['exonsDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_intronsDistance'] = np.where(df['intronsDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_promoterDistance'] = np.where(df['promoterDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_cgiDistance'] = np.where(df['cgiDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_ctcfDistance'] = np.where(df['ctcfDistance']<0, df['mixedReadCount'], 0)
    df['mixedReads_geneDistalRegulatoryModulesDistance'] = np.where(df['geneDistalRegulatoryModulesDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_vistaEnhancersDistance'] = np.where(df['vistaEnhancersDistance'] <0, df['mixedReadCount'], 0)
    df['mixedReads_firstExonDistance'] = np.where(df['firstExonDistance'] <0, df['mixedReadCount'], 0)
                                                                        
                                                                        
    df['fullMethReads_genesDistance'] = np.where(df['genesDistance']<0, df['methReadCount'], 0)
    df['fullMethReads_exonsDistance'] = np.where(df['exonsDistance']<0, df['methReadCount'], 0)
    df['fullMethReads_intronsDistance'] = np.where(df['intronsDistance']<0, df['methReadCount'], 0)
    df['fullMethReads_promoterDistance'] = np.where(df['promoterDistance']<0, df['methReadCount'], 0)
    df['fullMethReads_cgiDistance'] = np.where(df['cgiDistance']<0, df['methReadCount'], 0)
    df['fullMethReads_ctcfDistance'] = np.where(df['ctcfDistance']<0, df['methReadCount'], 0)
    df['fullMethReads_geneDistalRegulatoryModulesDistance'] = np.where(df['geneDistalRegulatoryModulesDistance'] <0, df['methReadCount'], 0)
    df['fullMethReads_vistaEnhancersDistance'] = np.where(df['vistaEnhancersDistance'] <0, df['methReadCount'], 0)
    df['fullMethReads_firstExonDistance'] = np.where(df['firstExonDistance'] <0, df['methReadCount'], 0)
                                                                                                            
    df = df.sum()
    df['filename'] = str(filename)
                                                                                                                
    df['PDR_total'] = df['mixedReadCount']/df['total_reads']
    df['PDR_GenesBody'] = df['mixedReads_genesDistance']/df['totreads_genesDistance']
    df['PDR_Exons'] = df['mixedReads_exonsDistance']/df['totreads_exonsDistance']
    df['PDR_Introns'] = df['mixedReads_intronsDistance']/df['totreads_intronsDistance']
    df['PDR_Promoters'] = df['mixedReads_promoterDistance']/df['totreads_promoterDistance']
    df['PDR_CGIslands'] = df['mixedReads_cgiDistance']/df['totreads_cgiDistance']
    df['PDR_CTCF'] = df['mixedReads_ctcfDistance']/df['totreads_ctcfDistance']
    df['PDR_Enhancer'] = df['mixedReads_geneDistalRegulatoryModulesDistance']/df['totreads_geneDistalRegulatoryModulesDistance']
                                                                                                                                                
    df['percent_totalMeth'] = df['methReadCount']/df['total_reads']
    df['totalMeth_GenesBody'] = df['fullMethReads_genesDistance']/df['totreads_genesDistance']
    df['totalMeth_Exons'] = df['fullMethReads_exonsDistance']/df['totreads_exonsDistance']
    df['totalMeth_Introns'] = df['fullMethReads_intronsDistance']/df['totreads_intronsDistance']
    df['totalMeth_Promoters'] = df['fullMethReads_promoterDistance']/df['totreads_promoterDistance']
    df['totalMeth_CGIslands'] = df['fullMethReads_cgiDistance']/df['totreads_cgiDistance']
    df['totalMeth_CTCF'] = df['fullMethReads_ctcfDistance']/df['totreads_ctcfDistance']
    df['totalMeth_Enhancer'] = df['fullMethReads_geneDistalRegulatoryModulesDistance']/df['totreads_geneDistalRegulatoryModulesDistance']
    newdf1 = newdf1.append(df, ignore_index=True)



# export as .csv
newdf1.to_csv('PDR_genomicRegions_RRBS_NormalBCD19pCD27mcell67_88.csv')

################################################
################################################
#
# NOTE: We need to recalculate for annotations
#
# Dataframe headers versus labels in Landau et al (2014)
#
# 'tssDistance'   // always seems to be a zero value---why?
# 'genesDistance' = 'Genes Body'
# 'exonsDistance' = 'Exons'
# 'intronsDistance' = 'Introns'
# 'promoterDistance' = 'Promoters'
# 'cgiDistance' = 'CG Islands'
# 'ctcfDistance' = 'CTCF binding site density'
# 'ctcfUpDistance'
# 'ctcfDownDistance'
# 'geneDistalRegulatoryModulesDistance' = 'Enhancer'
# 'vistaEnhancersDistance' = // ignore
# 'firstExonDistance'
#
###############
# QUESTIONS
###############
#
# Question (1) Calculating GCI shores and shelves is very tricky, as one must know the exact GCI boundaries
# e.g.
# if GCI distance is -1, this is included in GCIshore up [0 to 2000]
# if GCI distance is -2001, this is included in GCIshelf up [2000 to 4000]
#
# One cannot do this:
# df['GCIshoreUp'] = df['cgiDistance'] + 2000
# df['GCIshoreDown'] = df['cgiDistance'] - 2000
# df['GCIshelfUp'] = df['cgiDistance'] + 4000
# df['GCIshelfDown'] = df['cgiDistance'] - 4000
# as you are using 'cgiDistance' to be both the left boundary and the right boundary
#
# Question (2) How to calculate "Intergenic"?
#
# Question (3) What's up with 'tssDistance'?
#

