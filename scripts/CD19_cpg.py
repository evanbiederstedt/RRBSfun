import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import os
os.chdir('/Users/evanbiederstedt/Downloads/RRBS_data_files')



cd19cell = glob.glob("stacked_RRBS_NormalBCD19pcell*")

newdf1 = pd.DataFrame()
for filename in cd19cell:
    df = pd.read_csv(filename)
    df['filename'] = str(filename)
    df = df.sum()
    df["total_cpg_no_filter"] = df["avgReadCpGs"]
    df["filename"] = str(filename)
    df["filename"] = df["filename"][8:50]
    
    newdf1 = newdf1.append(df, ignore_index=True)

newdf1.to_csv("cd19_cpg1.csv")

cd19cell = glob.glob("stacked_RRBS_NormalBCD19pcell*")

newdf2 = pd.DataFrame()
for filename in cd19cell:
    df = pd.read_csv(filename)
    df['filename'] = str(filename)
    df = df[df['avgReadCpGs'] > 1]
    df = df.sum()
    df["total_cpg_gtrthan1"] = df["avgReadCpGs"]
    df["filename"] = str(filename)
    df["filename"] = df["filename"][8:50]
    newdf2 = newdf2.append(df, ignore_index=True)

newdf2.to_csv("cd19_cpg2.csv")

cd19cell = glob.glob("stacked_RRBS_NormalBCD19pcell*")


newdf3 = pd.DataFrame()
for filename in cd19cell:
    df = pd.read_csv(filename)
    df['filename'] = str(filename)
    df = df[df['avgReadCpGs'] >= 3.8]
    df = df.sum()
    df["total_cpg_gtrthan38"] = df["avgReadCpGs"]
    df["filename"] = str(filename)
    df["filename"] = df["filename"][8:50]
    newdf3 = newdf3.append(df, ignore_index=True)

newdf3.to_csv("cd19_cpg3.csv")


