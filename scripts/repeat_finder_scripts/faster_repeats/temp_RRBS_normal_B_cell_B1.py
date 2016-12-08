import glob
import pandas as pd
import numpy as np


df1 = pd.read_csv("repeats_hg19.csv")

RRBS_files = glob.glob("RRBS_normal_B_cell_B1*")

df_dict = {group : df for group, df in df1.groupby(by="chr")}


# In[11]:

from numpy import nan

def between_range(row, group_dict):
    # get sub dataframe from dictionary, if key not found return nan
    subset = group_dict.get(row['chr'], None)
    if subset is None:
        return ''

    # slice subset to get even smaller group
    subset = subset.loc[subset['start'] <= row['start'], :]

    # check if subset has any values
    if subset.empty:
        return ''

    # now slice the second time and return only a series of 'labels'
    subset = subset.loc[subset['end'] >= row['start'], 'repeat_class']

    # now you need to decide what to do if there are more than one labels, i have only taken here the first one
    # return the first label, if empty return nan
    if subset.empty:
        return ''
    else:
        return subset.iloc[0]


# In[12]:

from functools import partial
from time import time
between_range_partial = partial(between_range, group_dict=df_dict)

cols = ['chr', 'start']

for filename in RRBS_files:
    df2 = pd.read_csv(filename, sep="\t")
    labels = df2.loc[:, cols].apply(between_range_partial, axis=1)
    df2["repeat_class"] = labels
    df2.to_csv(filename, sep='\t', index=False)

