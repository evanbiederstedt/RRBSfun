# fifth method
import glob
import numpy as np
import pandas as pd
import os

df1 = pd.read_csv("repeats_hg19_begin.csv")

anno = pd.read_table("RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC.dan.anno")

chr_one = anno[anno.chr == "chr1"]

df = chr_one.merge(df1, how = "left")
df = df.loc[(df.start >= df.begin) & (df.start <= df.end), ["chr", "start", "repeat_class"]]
result = chr_one.merge(df, how = "left")

result.to_csv("repeat_chr1_test5.csv", index=False)

