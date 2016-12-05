import glob
import numpy as np
import pandas as pd
from numpy import nan
import os

os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/RRBS_anno_clean")

repeats = pd.read_csv("repeats_hg19.csv")
annofiles = glob.glob("RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.A*")

def between_range(row):
    subset = repeats.loc[(row["chr"] == repeats.chr) & (row.start >= repeats.start) & (row.start <= repeats.end), :]
    if subset.empty:
        return np.nan
    return subset.repeat_class

#newdf1 = pd.DataFrame()
for filename in annofiles:
    df = pd.read_table(filename)
    df["hg19_repeats"] = df.apply(between_range, axis = 1)
    df.to_csv(str("repeatregions_") + filename + ".csv", index=False)

