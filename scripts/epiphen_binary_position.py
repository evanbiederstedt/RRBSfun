import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
pd.set_option('display.max_columns', 50) # print all rows


import os
os.chdir("/Users/evanbiederstedt/Downloads/phylo_files")

all_files = glob.glob("RRBS*")

for filename in all_files:
    df= pd.read_table(filename)
    df['filename'] = str('filename')
    df = df.drop(df.ix[:,"bedCol.2": "hypoInHues64Distance"], axis=1) # remove columns
    df["coord"] = df["start"].apply(str)  # convert to string
    df["position"] =df[["chr", "coord"]].apply(lambda x: '_'.join(x), axis=1) # join these two columns to create position
    df = df.drop("bedCol", 1)
    df = df.drop("coord", 1)
    df = df.rename(columns = {'bedCol.1':'Methyl'})  # change into 'Methy
    df["Methyl"] = df["Methyl"].str.replace(r'[\'\"]', '')
    df["Methyl"] = df.Methyl.apply(eval)
    df.loc[df["Methyl"] >= 0.75, "Methyl"] = 1
    df.loc[df["Methyl"] <= 0.25, "Methyl"] = 0
    df.loc[(df["Methyl"] < 0.75) & (df["Methyl"] > 0.25), "Methyl"] = str("NaN")
    df = df[["position", "Methyl"]]
    df.to_csv(str("binary_position_") + str(filename) +str(".csv"))