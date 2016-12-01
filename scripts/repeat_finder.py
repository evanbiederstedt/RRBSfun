import glob
import numpy as np
import pandas as pd
import os
import dask.dataframe as dd

repeats = dd.read_csv("repeats_hg19.csv")

anno = dd.read_table("RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC.dan.anno")

df1 = dd.merge(anno, repeats, on="chr", how="outer", suffixes=("","_repeat"))
df1.to_csv("find_repeatsTESToutput.csv", index=False)
df1 = df1[(repeats.chr == row.chr) & (anno.start >= repeats.begin) & (anno.start <= repeats.end)]
df1 = dd.merge(anno, df1, on = ["chr"])
df1.to_csv("find_repeatsTEST2.csv", index=False).compute(num_workers=20)
