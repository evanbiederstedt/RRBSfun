import glob
import numpy as np
import pandas as pd
import os
import dask.dataframe as dd

repeats = dd.read_csv("RRBS_repeats.csv")

anno = dd.read_table("RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC.dan.anno")

df1 = dd.merge(anno, repeats, on="chr", how="outer", suffixes=("","_repeat")).compute(num_workers=20)
df1.to_csv("find_repeatsTEST1014.csv", index=False)
df1 = df1[(repeats.chr == row.chr) & (anno.start >= repeats.begin) & (anno.start <= repeats.end)]
df1 = dd.merge(anno, df1, on = ["chr"]).compute(num_workers=20)
df1.to_csv("find_repeats1014.csv", index=False)




mpirun -np 40 ./examl-AVX -t /gpfs/commons/home/biederstedte-934/evan_projects/standard-RAxML-master/RAxML_info.total_chrom1.tree -m GAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files/total_chrom1.phy -n T1001



mpirun -np 40 ./examl-AVX -t /gpfs/commons/home/biederstedte-934/evan_projects/standard-RAxML-master/RAxML_info.total_chrom2.tree -m GAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/ExaML-master/parser/chrom2.binary -n T2






#!/bin/bash

mpirun -np 10 ./examl-AVX -t /gpfs/commons/home/biederstedte-934/evan_projects/standard-RAxML-master/RAxML_parsimonyTree.total_chrom4.tree -m GAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/ExaML-master/parser/chrom4.binary -n T4


./examl-AVX -t /gpfs/commons/home/biederstedte-934/evan_projects/standard-RAxML-master/RAxML_parsimonyTree.total_chrom4.tree -m GAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/ExaML-master/parser/chrom4.binary -n T4
