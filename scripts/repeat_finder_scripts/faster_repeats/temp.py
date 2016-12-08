
import pandas as pd
import numpy as np


df1 = pd.read_csv("repeats_hg19.csv")
df2 = pd.read_csv("RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC.dan.anno", sep="\t")



df1.shape

df2.shape

df1.head()

df2.head()

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
labels = df2.loc[:, cols].apply(between_range_partial, axis=1)


df2["repeat_class"] = labels
filepath = "TESTRRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC_added_class.dan.anno.csv"
df2.to_csv(filepath, sep='\t', index=False)

