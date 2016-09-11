import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")


pcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27pcell*")
print(len(pcell))
cw154cell = glob.glob("binary_position_RRBS_cw154*")
print(len(cw154cell)

totalfiles = pcell + cw154cell
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

total_matrix.columns = ['RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACAACC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACCGCG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACTCAC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.AGGATG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ATAGCG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ATCGAC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CAAGAG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CATGAC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CCTTCG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CGGTAG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CTATTG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.CTCAGC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GCATTC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GCTGCC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GGCATC',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GTGAGG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.GTTGAG',
                        'RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.TAGCGG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACAACC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACCGCG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACGTGG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ACTCAC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.AGGATG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ATAGCG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.ATCGAC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CAAGAG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CATGAC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CCTTCG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CGGTAG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CTATTG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.CTCAGC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GACACG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GCATTC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GCTGCC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GGCATC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GTGAGG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.GTTGAG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.TAGCGG',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.TATCTC',
                        'RRBS_NormalBCD19pCD27pcell23_44_CTCTCTAC.TCTCTG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.ACCGCG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.ACTCAC',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.ATAGCG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.CAAGAG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.CCTTCG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.CTATTG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.GACACG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.GTGAGG',
                        'RRBS_NormalBCD19pCD27pcell45_66_CAGAGAGG.TAGCGG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACAACC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACCGCG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACGTGG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ACTCAC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.AGGATG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ATAGCG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.ATCGAC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CATGAC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CCTTCG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CGGTAG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CTATTG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.CTCAGC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GACACG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GCATTC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GCTGCC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GGCATC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GTGAGG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.GTTGAG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.TAGCGG',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.TATCTC',
                        'RRBS_NormalBCD19pCD27pcell67_88_GCTACGCT.TCTCTG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACAACC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACCGCG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACGTGG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACTCAC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.AGGATG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ATAGCG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ATCGAC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CAAGAG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CATGAC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CCTTCG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CGGTAG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CTCAGC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GACACG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GCATTC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GCTGCC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GGCATC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GTGAGG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.TAGCGG',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.TATCTC',
                        'RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.TCTCTG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.ACAACC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.ACCGCG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.ACGTGG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.ACTCAC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.AGGATG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.ATAGCG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.ATCGAC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.CATGAC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.CCTTCG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.CGGTAG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.CTATTG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.CTCAGC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.GACACG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.GCATTC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.GCTGCC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.GGCATC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.GTGAGG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.GTTGAG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.TAGCGG',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.TATCTC',
                        'RRBS_cw154_Tris_protease_CTCTCTAC.TCTCTG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACAACC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACCGCG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACGTGG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACTCAC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.AGGATG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.ATAGCG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.ATCGAC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.CATGAC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.CCTTCG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.CGGTAG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.CTATTG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.CTCAGC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.GACACG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.GCATTC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.GCTGCC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.GGCATC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.GTGAGG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.GTTGAG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.TAGCGG',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.TATCTC',
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.TCTCTG']

print(total_matrix.shape)



total_matrix = total_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))



total_matrix = total_matrix.astype(str).apply(''.join)




tott = pd.Series(total_matrix.index.astype(str).str.cat(total_matrix.astype(str),'    '))



tott.to_csv("normal_cll_3.phy", header=None, index=None)

print(tott.shape)

             

