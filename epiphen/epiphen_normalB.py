import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
pd.set_option('display.max_columns', 50) # print all rows

import os
os.chdir('/Users/evanbiederstedt/Downloads/required_binary_phylo_files')

import statsmodels.api as sm

normalB_files = glob.glob("binary_position_RRBS_normal_B_cell*")

df_list = []
for file in normalB_files:
    df = pd.read_csv(file)
    df = df.drop("Unnamed: 0", axis=1)
    df_list.append(df)

normalB_matrix = pd.concat([df.set_index("position") for df in df_list], axis=1).reset_index().astype(object)

normalB_matrix = normalB_matrix.drop("index", axis=1)

normalB_matrix.columns = ["RRBS_normal_B_cell_A1_24_TAAGGCGA.ACAACC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.ACCGCG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.ACGTGG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.AGGATG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.ATAGCG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.ATCGAC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.CAAGAG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.CATGAC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.CGGTAG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.CTATTG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.CTCAGC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.GACACG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.GCTGCC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.GGCATC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.GTGAGG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.GTTGAG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.TAGCGG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.TATCTC",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.TCTCTG",
                          "RRBS_normal_B_cell_A1_24_TAAGGCGA.TGACAG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.ACAACC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.ACCGCG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.ACTCAC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.ATAGCG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.CAAGAG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.CATGAC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.CCTTCG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.CGGTAG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.CTATTG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.CTCAGC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.GACACG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.GCATTC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.GGCATC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.GTGAGG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.GTTGAG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.TAGCGG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.TATCTC",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.TCTCTG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.TGACAG",
                          "RRBS_normal_B_cell_B1_24_CGTACTAG.TGCTGC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.ACAACC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.ACCGCG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.ACGTGG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.ACTCAC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.AGGATG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.ATAGCG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.ATCGAC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.CAAGAG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.CATGAC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.CGGTAG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.CTATTG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.GACACG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.GCATTC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.GCTGCC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.GGCATC",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.GTGAGG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.GTTGAG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.TAGCGG",
                          "RRBS_normal_B_cell_C1_24_AGGCAGAA.TATCTC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.ACAACC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.ACCGCG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.ACGTGG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.ACTCAC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.AGGATG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.ATCGAC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.CAAGAG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.CATGAC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.CCTTCG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.CGGTAG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.CTATTG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.CTCAGC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.GACACG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.GCATTC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.GCTGCC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.GGCATC",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.GTTGAG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.TAGCGG",
                          "RRBS_normal_B_cell_D1_24_TCCTGAGC.TATCTC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.ACAACC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.ACCGCG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.ACGTGG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.ACTCAC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.AGGATG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.ATAGCG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.ATCGAC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.CAAGAG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.CATGAC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.CGGTAG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.CTATTG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.CTCAGC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.GACACG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.GCATTC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.GCTGCC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.GGCATC",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.GTGAGG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.TAGCGG",
                          "RRBS_normal_B_cell_G1_22_GGACTCCT.TATCTC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.ACCGCG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.ACGTGG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.ACTCAC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.AGGATG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.ATCGAC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.CAAGAG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.CATGAC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.CCTTCG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.CTATTG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.CTCAGC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.GCATTC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.GCTGCC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.GGCATC",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.GTGAGG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.GTTGAG",
                          "RRBS_normal_B_cell_H1_22_TAGGCATG.TCTCTG"]

normalB_matrix = normalB_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))

normalB_matrix = normalB_matrix.astype(str).apply(''.join)

tot_normB = pd.Series(normalB_matrix.index.astype(str).str.cat(normalB_matrix.astype(str),'    '))

tot_normB.to_csv("RRBS_NormalB_total.phy", header=None, index=None)








