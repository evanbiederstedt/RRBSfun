import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
pd.set_option('display.max_columns', 50) # print all rows

# proof of concept

import os
os.chdir("/Users/evanbiederstedt/Downloads/phylo_files")

df = pd.read_table("RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GTTGAG.bed.anno")

df = df.drop(df.ix[:,'bedCol.2': 'hypoInHues64Distance'], axis=1)
df = df.drop("bedCol", 1)
df = df.rename(columns = {'bedCol.1':"Methyl"})

df["Methyl"] = df["Methyl"].str.replace(r'[\'\"]', '')
df["Methyl"] = df.Meth.apply(eval)
df.loc[df["Methyl"] >= 0.75, "Methyl"] = 1
df.loc[df["Methyl"] <= 0.25, "Methyl"] = 0
df.loc[(df["Methyl"] < 0.75) & (df["Methyl"] > 0.25), "Methyl"] = str('NaN')

df.shape



