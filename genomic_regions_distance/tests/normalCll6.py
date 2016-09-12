import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")

mcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27mcell*")
print(len(mcell))
cw154cell = glob.glob("binary_position_RRBS_cw154*")
print(len(cw154cell))

totalfiles = mcell + cw154cell
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

total_matrix.columns = ['RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ACAACC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ACGTGG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ACTCAC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ATAGCG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ATCGAC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CAAGAG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CATGAC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CCTTCG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CTATTG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CTCAGC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GACACG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GCATTC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GCTGCC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GGCATC',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GTGAGG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GTTGAG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.TAGCGG',
                              'RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.TATCTC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACAACC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACCGCG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACGTGG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACTCAC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.AGGATG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ATAGCG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ATCGAC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CAAGAG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CATGAC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CCTTCG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CGGTAG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CTATTG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CTCAGC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GACACG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GCATTC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GTGAGG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GTTGAG',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.TATCTC',
                              'RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.TCTCTG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ACAACC',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ACGTGG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ACTCAC',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.AGGATG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ATAGCG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ATCGAC',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CAAGAG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CATGAC',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CCTTCG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CGGTAG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CTATTG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CTCAGC',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.GACACG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.GTGAGG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.TAGCGG',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.TATCTC',
                              'RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.TCTCTG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACAACC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACCGCG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACGTGG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACTCAC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.AGGATG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ATAGCG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ATCGAC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CAAGAG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CATGAC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CCTTCG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CGGTAG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CTATTG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CTCAGC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GACACG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GCATTC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GGCATC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GTGAGG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GTTGAG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.TAGCGG',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.TATCTC',
                              'RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.TCTCTG',
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



tott.to_csv("normal_cll_6.phy", header=None, index=None)

print(tott.shape)

             

