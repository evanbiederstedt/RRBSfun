import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
pd.set_option('display.max_columns', 50) # print all rows


import os
os.chdir("/Users/evanbiederstedt/Downloads/phylo_files")

cw154 = glob.glob("RRBS_cw154*")

for filename in cw154:
    df= pd.read_table(filename)
    df['filename'] = str('filename')
    df = df.drop(df.ix[:,'bedCol.2': 'hypoInHues64Distance'], axis=1)
    df = df.drop("bedCol", 1)
    df = df.rename(columns = {'bedCol.1':'Methyl'})
    df["Methyl"] = df["Methyl"].str.replace(r'[\'\"]', '')
    df["Meth"] = df.Meth.apply(eval)
    df.loc[df['Meth'] >= 0.75, 'Meth'] = 1
    df.loc[df['Meth'] <= 0.25, 'Meth'] = 0
    df.loc[(df["Meth"] < 0.75) & (df["Meth"] > 0.25), 'Meth'] = str('NaN')
    df.to_csv(str("phyl_") + str(filename) +str(".csv"))