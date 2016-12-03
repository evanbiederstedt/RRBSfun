# fourth method
import glob
import numpy as np
import pandas as pd
import os

df1 = pd.read_csv("repeats_hg19.csv")

anno = pd.read_table("RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC.dan.anno")

def match_labels(row):
    current_df = df1[(df1["chr"]==row["chr"]) & (df1["start"]<=row["start"]) & (df1["end"]>=row["start"])]
    try:
        row["repeat_class"] = current_df["repeat_class"].iloc[0]
    except:
        row["repeat_class"] = np.NaN
    return row

result = anno.apply(lambda x:match_labels(x), axis=1)

result.to_csv("results_repeats4.csv", index=False)
