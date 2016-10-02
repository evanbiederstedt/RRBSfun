import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/projects_evan/correct_phylo_files")


cw154cell = glob.glob("binary_position_RRBS_cw154_Tris_protease*")
print(len(cw154cell))
cd19cell = glob.glob("binary_position_RRBS_NormalBCD19pcell*")
print(len(cd19cell))

totalfiles = cw154cell + cd19cell
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

total_matrix.columns = ["RRBS_cw154_Tris_protease_CTCTCTAC.ACAACC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACCGCG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACGTGG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACTCAC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.AGGATG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ATAGCG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ATCGAC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CATGAC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CCTTCG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CGGTAG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CTATTG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CTCAGC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GACACG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GCATTC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GCTGCC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GGCATC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GTGAGG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GTTGAG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.TAGCGG"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.TATCTC"
                        "RRBS_cw154_Tris_protease_CTCTCTAC.TCTCTG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACAACC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACCGCG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACGTGG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACTCAC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.AGGATG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ATAGCG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ATCGAC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CATGAC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CCTTCG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CGGTAG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CTATTG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CTCAGC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GACACG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GCATTC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GCTGCC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GGCATC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GTGAGG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GTTGAG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.TAGCGG"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.TATCTC"
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.TCTCTG"
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.ACAACC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.ACCGCG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.ACGTGG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.ACTCAC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.AGGATG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.ATAGCG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.ATCGAC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.CAAGAG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.CATGAC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.CCTTCG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.CGGTAG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.CTATTG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.CTCAGC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GACACG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GCATTC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GCTGCC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GGCATC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GTGAGG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GTTGAG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TAGCGG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TATCTC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TCTCTG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TGCTGC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACAACC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACCGCG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACGTGG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACTCAC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.AGGATG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ATAGCG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ATCGAC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.CAAGAG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.CATGAC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.CCTTCG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.CGGTAG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.CTATTG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.CTCAGC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.GACACG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.GCATTC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.GCTGCC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.GGCATC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.GTGAGG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.GTTGAG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TAGCGG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TATCTC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TCTCTG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TGCTGC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.ACAACC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.ACCGCG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.ACGTGG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.ACTCAC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.AGGATG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.ATAGCG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.ATCGAC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.CAAGAG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.CATGAC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.CCTTCG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.CGGTAG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.CTATTG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.CTCAGC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.GACACG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.GCATTC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.GCTGCC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.GGCATC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.GTGAGG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.GTTGAG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.TAGCGG",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.TATCTC",
                        "RRBS_NormalBCD19pcell45_66_AGGCAGAA.TCTCTG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.ACAACC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.ACCGCG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.ACGTGG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.ACTCAC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.AGGATG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.ATAGCG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.ATCGAC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.CAAGAG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.CATGAC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.CCTTCG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.CGGTAG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.CTATTG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.CTCAGC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GACACG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GCATTC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GCTGCC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GGCATC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GTGAGG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GTTGAG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TAGCGG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TATCTC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TCTCTG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TGCTGC"]

print(total_matrix.shape)



total_matrix = total_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))



total_matrix = total_matrix.astype(str).apply(''.join)




tott = pd.Series(total_matrix.index.astype(str).str.cat(total_matrix.astype(str),'    '))



tott.to_csv("normal_cll_2.phy", header=None, index=None)

print(tott.shape)

            
