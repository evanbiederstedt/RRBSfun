import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")


normalb = glob.glob("binary_position_RRBS_normal_B_cell*")
print(len(normalb))
cw154cell = glob.glob("binary_position_RRBS_cw154*")
print(len(cw154cell))


totalfiles = normalb + cw154cell
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

total_matrix.columns = ['RRBS_normal_B_cell_A1_24_TAAGGCGA.ACAACC',
                        'RRBS_normal_B_cell_A1_24_TAAGGCGA.ACCGCG',
                        'RRBS_normal_B_cell_A1_24_TAAGGCGA.ACGTGG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.AGGATG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.ATAGCG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.ATCGAC',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.CAAGAG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.CATGAC',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.CGGTAG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.CTATTG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.CTCAGC',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.GACACG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.GCTGCC',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.GGCATC',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.GTGAGG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.GTTGAG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.TAGCGG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.TATCTC',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.TCTCTG',
                              'RRBS_normal_B_cell_A1_24_TAAGGCGA.TGACAG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.ACAACC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.ACCGCG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.ACTCAC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.ATAGCG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.CAAGAG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.CATGAC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.CCTTCG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.CGGTAG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.CTATTG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.CTCAGC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.GACACG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.GCATTC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.GGCATC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.GTGAGG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.GTTGAG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.TAGCGG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.TATCTC',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.TCTCTG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.TGACAG',
                              'RRBS_normal_B_cell_B1_24_CGTACTAG.TGCTGC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.ACAACC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.ACCGCG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.ACGTGG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.ACTCAC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.AGGATG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.ATAGCG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.ATCGAC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.CAAGAG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.CATGAC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.CGGTAG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.CTATTG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.GACACG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.GCATTC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.GCTGCC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.GGCATC',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.GTGAGG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.GTTGAG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.TAGCGG',
                              'RRBS_normal_B_cell_C1_24_AGGCAGAA.TATCTC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.ACAACC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.ACCGCG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.ACGTGG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.ACTCAC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.AGGATG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.ATCGAC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.CAAGAG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.CATGAC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.CCTTCG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.CGGTAG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.CTATTG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.CTCAGC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.GACACG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.GCATTC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.GCTGCC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.GGCATC',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.GTTGAG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.TAGCGG',
                              'RRBS_normal_B_cell_D1_24_TCCTGAGC.TATCTC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.ACAACC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.ACCGCG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.ACGTGG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.ACTCAC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.AGGATG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.ATAGCG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.ATCGAC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.CAAGAG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.CATGAC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.CGGTAG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.CTATTG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.CTCAGC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.GACACG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.GCATTC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.GCTGCC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.GGCATC',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.GTGAGG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.TAGCGG',
                              'RRBS_normal_B_cell_G1_22_GGACTCCT.TATCTC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.ACCGCG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.ACGTGG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.ACTCAC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.AGGATG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.ATCGAC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.CAAGAG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.CATGAC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.CCTTCG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.CTATTG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.CTCAGC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.GCATTC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.GCTGCC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.GGCATC',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.GTGAGG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.GTTGAG',
                              'RRBS_normal_B_cell_H1_22_TAGGCATG.TCTCTG',
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



tott.to_csv("normal_cll_8.phy", header=None, index=None)

print(tott.shape)

             

