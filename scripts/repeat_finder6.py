# third approach
import glob
import numpy as np
import pandas as pd
import os
from functools import partial
from numpy import nan

df1 = pd.read_csv("repeats_hg19.csv")

anno = pd.read_table("RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC.dan.anno")

# first split anno into subgroups with `groupby`
df_dict = {group : df for group, df in anno.groupby(by="chr")}

#
# outputs keys
# ['chr12', 'chr14', 'chr5', 'chr21', 'chr20', 'chr17', 'chr2',
# 'chr22', 'chr10', 'chr11', 'chr13', 'chr6', 'chr19', 'chr16',
# 'chr4', 'chr15', 'chr1', 'chr8', 'chr7', 'chr9', 'chr3', 'chr18']
#

def between_range(row, group_dict):
    # searching by each "chr"
    #
    # get sub dataframe from dictionary, if key not found return np.NaN
    subset = group_dict.get(row["chr"], None)
    if subset is None:
        return ''

    # slice subset to get even smaller group
    subset = subset.loc[subset["start"] <= row["start"], :]
    # check if subset has any values
    if subset.empty:
        return ''

    # now slice the second time and return only a series of 'labels'
    subset = subset.loc[subset["end"] >= row["start"], "repeat_class"]

    # return the first "repeat_class" category; otherwise
    if subset.empty:
        return ''
    else:
        return subset.iloc[0]

# now use the magic of functools to get a pandas Series
#   that you can add as a new column for the second df "anno"

between_range_partial = partial(between_range, group_dict=df_dict)

cols = ["chr", "start"]

repeats_found = anno.loc[:, cols].apply(between_range_partial, axis=1)
anno["hg19_repeats"] = repeats_found

anno.to_csv("repeats_try6.csv", index=False)




