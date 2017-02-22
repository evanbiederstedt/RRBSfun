#!/bin/python3

import csv
import pandas as pd
from collections import Counter

results = {}
with open("total_cells.phy") as f:
    next(f)  # skip first line
    next(f)
    for line in f:
        title, code = line.strip().split()
        results[title] = Counter(line)

## remove superfluous keys in dicitonary

keep_keys = ["?", "1", "0"]

def filter_dict_results(items):
    tmp_dict = [{k : v} for k,v in items if k in keep_keys]
    ret = {}
    for elem in tmp_dict:
        ret.update(elem)
    return [ret]

corrected_dictionary1 = [ (k, res_dict) for (k, v) in results.items() for res_dict in filter_dict_results(v.items())]

print(corrected_dictionary1)

df = pd.DataFrame.from_dict(dict(corrected_dictionary1), orient="index").rename_axis("cells").reset_index()

df.to_csv("total_cells_matrix_counts.csv", index=False)

