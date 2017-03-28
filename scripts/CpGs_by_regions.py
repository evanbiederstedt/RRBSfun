#!/usr/bin/python3

import glob
import pandas as pd
import numpy as np
import os


os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/RRBS_anno_clean4")


all_files = glob.glob("*dan.anno")

newdf1 = pd.DataFrame()
for filename in all_files:
    df = pd.read_table(filename)
    df["filename"] = str(filename)
    ### drop unnecessary columns for this task
    df = df.drop(["chr", "start", "strand", "avgWeightedEnt", "CpGEntropy", "tss", "genes", "exons", "introns", "promoter",
                  "cgi", "geneDensity", "ctcf", "ctcfUpstream", "ctcfDownstream", "ctcfDensity", "geneDistalRegulatoryModules",
                  "vistaEnhancers","3PrimeUTR",  "5PrimeUTR", "firstExon", "geneDistalRegulatoryModulesK562","hypoInHues64"], axis=1)
     
    df["avgReadCpGs"] = np.round(df["avgReadCpGs"])
    ### create unique Stack ID for each set of rows with the shared CpGs
    num_bins2 = np.ceil(df["avgReadCpGs"].max()/1.25)
                  
    df["avgReadCpGs_binned"] = pd.cut(df["avgReadCpGs"], num_bins2, labels=False)
                  
    df["read_stack_ID"] = (df.avgReadCpGs_binned.shift(1) != df.avgReadCpGs_binned).astype(int).cumsum()

    df = df.drop_duplicates("read_stack_ID")

    intergenic_set = ( (df["genesDistance"]>0) & (df["exonsDistance"]>0) & (df["intronsDistance"]>0) & (df["promoterDistance"]>0) & (df["cgiDistance"]>0) & (df["ctcfDistance"]>0) &
                  (df["ctcfUpDistance"]>0) & (df["ctcfDownDistance"]>0) & (df["geneDistalRegulatoryModulesDistance"]>0) & (df["vistaEnhancersDistance"]>0) &
                  (df["3PrimeUTRDistance"]>0) & (df["5PrimeUTRDistance"]>0) & (df["firstExonDistance"]>0) & (df["geneDistalRegulatoryModulesK562Distance"]>0) & (df["hypoInHues64Distance"]>0) )

    ##
    ## shore (+1 to +2000....on either side)
    ## shelf (+2000 to +4000....on either side)
    ##

    df["avgReadCpGs_genesDistance"] = np.where(df["genesDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_exonsDistance"] = np.where(df["exonsDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_intronsDistance"] = np.where(df["intronsDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_promoterDistance"] = np.where(df["promoterDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_cgiDistance"] = np.where(df["cgiDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_ctcfDistance"] = np.where(df["ctcfDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_ctcfUpDistance"] = np.where(df["ctcfUpDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_ctcfDownDistance"] = np.where(df["ctcfDownDistance"]<0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_geneDistalRegulatoryModulesDistance"] = np.where(df["geneDistalRegulatoryModulesDistance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_vistaEnhancersDistance"] = np.where(df["vistaEnhancersDistance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_3PrimeUTRDistance"] = np.where(df["3PrimeUTRDistance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_5PrimeUTRDistance"] = np.where(df["5PrimeUTRDistance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_firstExonDistance"] = np.where(df["firstExonDistance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_geneDistalRegulatoryModulesK562Distance"] = np.where(df["geneDistalRegulatoryModulesK562Distance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_hypoInHues64Distance"] = np.where(df["hypoInHues64Distance"] <0, df["avgReadCpGs"], 0)
    df["avgReadCpGs_intergenic"] = np.where(intergenic_set, df["avgReadCpGs"], 0)
    df["avgReadCpGs_shore"] = np.where((df["cgiDistance"].between(0, 2000, inclusive=False)), df["avgReadCpGs"], 0)
    df["avgReadCpGs_shelf"] = np.where(((df["cgiDistance"] >= 2000) & (df["cgiDistance"] < 4000)), df["avgReadCpGs"], 0)

    ### take sum
    df["total_avgReadCpGs_genesDistance"]  = df["avgReadCpGs_genesDistance"].sum()
    df["total_avgReadCpGs_exonsDistance"]  = df["avgReadCpGs_exonsDistance"].sum()
    df["total_avgReadCpGs_intronsDistance"] = df["avgReadCpGs_intronsDistance"].sum()
    df["total_avgReadCpGs_promoterDistance"] = df["avgReadCpGs_promoterDistance"].sum()
    df["total_avgReadCpGs_cgiDistance"] = df["avgReadCpGs_cgiDistance"].sum()
    df["total_avgReadCpGs_ctcfDistance"] = df["avgReadCpGs_ctcfDistance"].sum()
    df["total_avgReadCpGs_ctcfUpDistance"] = df["avgReadCpGs_ctcfUpDistance"].sum()
    df["total_avgReadCpGs_ctcfDownDistance"] = df["avgReadCpGs_ctcfDownDistance"].sum()
    df["total_avgReadCpGs_geneDistalRegulatoryModulesDistance"] = df["avgReadCpGs_geneDistalRegulatoryModulesDistance"].sum()
    df["total_avgReadCpGs_vistaEnhancersDistance"] = df["avgReadCpGs_vistaEnhancersDistance"].sum()
    df["total_avgReadCpGs_3PrimeUTRDistance"] = df["avgReadCpGs_3PrimeUTRDistance"].sum()
    df["total_avgReadCpGs_5PrimeUTRDistance"] = df["avgReadCpGs_5PrimeUTRDistance"].sum()
    df["total_avgReadCpGs_firstExonDistance"] = df["avgReadCpGs_firstExonDistance"].sum()
    df["total_avgReadCpGs_geneDistalRegulatoryModulesK562Distance"] = df["avgReadCpGs_geneDistalRegulatoryModulesK562Distance"].sum()
    df["total_avgReadCpGs_hypoInHues64Distance"] = df["avgReadCpGs_hypoInHues64Distance"].sum()
    df["total_avgReadCpGs_intergenic"] = df["avgReadCpGs_intergenic"].sum()
    df["total_avgReadCpGs_shore"] = df["avgReadCpGs_shore"].sum()
    df["total_avgReadCpGs_shelf"] = df["avgReadCpGs_shelf"].sum()
                                                                             
    df["filename"] = str(filename)
    
    newdf1 = newdf1.append(df, ignore_index=True)

newdf1 = newdf1[["filename", "total_avgReadCpGs_genesDistance", "total_avgReadCpGs_exonsDistance", "total_avgReadCpGs_intronsDistance", "total_avgReadCpGs_promoterDistance", "total_avgReadCpGs_cgiDistance", "total_avgReadCpGs_ctcfDistance", "total_avgReadCpGs_ctcfUpDistance", "total_avgReadCpGs_ctcfDownDistance", "total_avgReadCpGs_geneDistalRegulatoryModulesDistance", "total_avgReadCpGs_vistaEnhancersDistance", "total_avgReadCpGs_3PrimeUTRDistance", "total_avgReadCpGs_5PrimeUTRDistance", "total_avgReadCpGs_firstExonDistance", "total_avgReadCpGs_geneDistalRegulatoryModulesK562Distance", "total_avgReadCpGs_hypoInHues64Distance", "total_avgReadCpGs_intergenic", "total_avgReadCpGs_shore", "total_avgReadCpGs_shelf"]]


newdf1.to_csv("CpGs_by_regions.csv", index=False)








