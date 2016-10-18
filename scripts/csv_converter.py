import glob
import numpy as np
import pandas as pd
import os


os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/RRBS_data_csv")

for file in glob.glob("RRBS*"):
    df = pd.read_table(file)
    df.to_csv(str(file) +".csv", index=False)

