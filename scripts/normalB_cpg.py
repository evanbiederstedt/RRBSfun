import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')



normal_B = glob.glob("stacked_RRBS_normal_B*")

newdf1 = pd.DataFrame()
for filename in normal_B[:2]:
    df = pd.read_csv(filename)
    df['filename'] = str(filename)
    df = df.sum()
    df["total_cpg_no_filter"] = df["avgReadCpGs"]
    df["filename"] = str(filename)
    df["filename"] = df["filename"][8:48]
    
    newdf1 = newdf1.append(df, ignore_index=True)

newdf1.to_csv("cpg_normalB1.csv")

normal_B = glob.glob("stacked_RRBS_normal_B*")

newdf2 = pd.DataFrame()
for filename in normal_B:
    df = pd.read_csv(filename)
    df['filename'] = str(filename)
    df = df[df['avgReadCpGs'] > 1]
    df = df.sum()
    df["total_cpg_gtrthan1"] = df["avgReadCpGs"]
    df["filename"] = str(filename)
    df["filename"] = df["filename"][8:48]
    newdf2 = newdf2.append(df, ignore_index=True)

newdf2.to_csv("cpg_normalB2.csv")



normal_B = glob.glob("stacked_RRBS_normal_B*")


newdf3 = pd.DataFrame()
for filename in normal_B:
    df = pd.read_csv(filename)
    df['filename'] = str(filename)
    df = df[df['avgReadCpGs'] >= 3.8]
    df = df.sum()
    df["total_cpg_gtrthan38"] = df["avgReadCpGs"]
    df["filename"] = str(filename)
    df["filename"] = df["filename"][8:48]
    newdf3 = newdf3.append(df, ignore_index=True)

newdf3.to_csv("cpg_normalB3.csv")



