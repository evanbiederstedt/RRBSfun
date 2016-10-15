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
