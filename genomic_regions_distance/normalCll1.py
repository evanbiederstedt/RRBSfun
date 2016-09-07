import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/personal_projects/correct_phylo_files")


pcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27pcell*")
print(len(pcell))
tritocell = glob.glob("binary_position_RRBS_trito_pool*")
print(len(tritocell))

totalfiles = pcell + tritocell
print(len(totalfiles))

df_list = []
for file in totalfiles:
    df = pd.read_csv(file)
    df = df.drop("Unnamed: 0", axis=1)
    df_list.append(df)


print(len(df_list))


total_matrix = pd.concat([df.set_index("position") for df in df_list], axis=1).reset_index().astype(object)

total_matrix = total_matrix.drop("index", axis=1)

len(total_matrix.columns)

total_matrix.columns = ["RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACAACC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACCGCG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACTCAC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.AGGATG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ATAGCG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ATCGAC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CAAGAG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CCTTCG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CGGTAG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CTATTG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CTCAGC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GCATTC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GCTGCC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GGCATC",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GTGAGG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GTTGAG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.TAGCGG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACAACC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACCGCG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACGTGG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACTCAC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.AGGATG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ATAGCG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ATCGAC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CAAGAG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CATGAC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CCTTCG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CGGTAG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CTATTG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CTCAGC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GACACG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GCATTC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GCTGCC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GGCATC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GTGAGG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GTTGAG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.TAGCGG",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.TATCTC",
                        "RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.TCTCTG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.ACCGCG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.ACTCAC",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.ATAGCG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.CAAGAG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.CCTTCG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.CTATTG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.GACACG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.GTGAGG",
                        "RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.TAGCGG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACAACC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACCGCG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACGTGG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACTCAC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.AGGATG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ATAGCG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ATCGAC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CATGAC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CCTTCG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CGGTAG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CTATTG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CTCAGC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GACACG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GCATTC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GCTGCC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GGCATC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GTGAGG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GTTGAG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.TAGCGG",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.TATCTC",
                        "RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.TCTCTG",
                        "RRBS_trito_pool_1_TAAGGCGA.ACAACC",
                        "RRBS_trito_pool_1_TAAGGCGA.ACGTGG",
                        "RRBS_trito_pool_1_TAAGGCGA.ACTCAC",
                        "RRBS_trito_pool_1_TAAGGCGA.ATAGCG",
                        "RRBS_trito_pool_1_TAAGGCGA.ATCGAC",
                        "RRBS_trito_pool_1_TAAGGCGA.CAAGAG",
                        "RRBS_trito_pool_1_TAAGGCGA.CATGAC",
                        "RRBS_trito_pool_1_TAAGGCGA.CCTTCG",
                        "RRBS_trito_pool_1_TAAGGCGA.CGGTAG",
                        "RRBS_trito_pool_1_TAAGGCGA.CTATTG",
                        "RRBS_trito_pool_1_TAAGGCGA.GACACG",
                        "RRBS_trito_pool_1_TAAGGCGA.GCATTC",
                        "RRBS_trito_pool_1_TAAGGCGA.GCTGCC",
                        "RRBS_trito_pool_1_TAAGGCGA.GGCATC",
                        "RRBS_trito_pool_1_TAAGGCGA.GTGAGG",
                        "RRBS_trito_pool_1_TAAGGCGA.GTTGAG",
                        "RRBS_trito_pool_1_TAAGGCGA.TAGCGG",
                        "RRBS_trito_pool_1_TAAGGCGA.TATCTC",
                        "RRBS_trito_pool_1_TAAGGCGA.TCTCTG",
                        "RRBS_trito_pool_1_TAAGGCGA.TGACAG",
                        "RRBS_trito_pool_1_TAAGGCGA.TGCTGC",
                        "RRBS_trito_pool_2_CGTACTAG.ACAACC",
                        "RRBS_trito_pool_2_CGTACTAG.ACGTGG",
                        "RRBS_trito_pool_2_CGTACTAG.ACTCAC",
                        "RRBS_trito_pool_2_CGTACTAG.AGGATG",
                        "RRBS_trito_pool_2_CGTACTAG.ATAGCG",
                        "RRBS_trito_pool_2_CGTACTAG.ATCGAC",
                        "RRBS_trito_pool_2_CGTACTAG.CAAGAG",
                        "RRBS_trito_pool_2_CGTACTAG.CATGAC",
                        "RRBS_trito_pool_2_CGTACTAG.CCTTCG",
                        "RRBS_trito_pool_2_CGTACTAG.CGGTAG",
                        "RRBS_trito_pool_2_CGTACTAG.CTATTG",
                        "RRBS_trito_pool_2_CGTACTAG.GACACG",
                        "RRBS_trito_pool_2_CGTACTAG.GCATTC",
                        "RRBS_trito_pool_2_CGTACTAG.GCTGCC",
                        "RRBS_trito_pool_2_CGTACTAG.GGCATC",
                        "RRBS_trito_pool_2_CGTACTAG.GTGAGG",
                        "RRBS_trito_pool_2_CGTACTAG.GTTGAG",
                        "RRBS_trito_pool_2_CGTACTAG.TAGCGG",
                        "RRBS_trito_pool_2_CGTACTAG.TATCTC",
                        "RRBS_trito_pool_2_CGTACTAG.TCTCTG",
                        "RRBS_trito_pool_2_CGTACTAG.TGACAG"]

print(total_matrix.shape)



total_matrix = total_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))



total_matrix = total_matrix.astype(str).apply(''.join)




tott = pd.Series(total_matrix.index.astype(str).str.cat(total_matrix.astype(str),'    '))



tott.to_csv("normal_cll_1.phy", header=None, index=None)

print(tott.shape)

             

