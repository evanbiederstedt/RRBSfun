import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 50) # print all rows


import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")


normalB = glob.glob("binary_position_RRBS_normal_B_cell*")
mcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27mcell*")
pcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27pcell*")
cd19cell = glob.glob("binary_position_RRBS_NormalBCD19pcell*")
cw154 = glob.glob("binary_position_RRBS_cw154*")
trito = glob.glob("binary_position_RRBS_trito_pool*")


print(len(normalB))
print(len(mcell))
print(len(pcell))
print(len(cd19cell))
print(len(cw154))
print(len(trito))

totalfiles = normalB + mcell + pcell + cd19cell + cw154 + trito
print(len(totalfiles))

df_list = []
for file in totalfiles:
    df = pd.read_csv(file)
    df = df.drop("Unnamed: 0", axis=1)
    df["chromosome"] = df["position"].map(lambda x: str(x)[:5])
    df = df[df["chromosome"] == "chr17"]
    df = df.drop("chromosome", axis=1)
    df_list.append(df)


print(len(df_list))


total_matrix = pd.concat([df.set_index("position") for df in df_list], axis=1).reset_index().astype(object)

total_matrix = total_matrix.drop("index", axis=1)

len(total_matrix.columns)


total_matrix.columns = ["RRBS_normal_B_cell_A1_24_TAAGGCGA.ACAACC",
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
                        "RRBS_normal_B_cell_H1_22_TAGGCATG.TCTCTG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ACCGCG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ACGTGG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ACTCAC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ATAGCG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.ATCGAC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CAAGAG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CATGAC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CCTTCG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CTATTG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.CTCAGC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GACACG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GCATTC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GCTGCC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GGCATC",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GTGAGG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.GTTGAG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.TAGCGG",
                        "RRBS_NormalBCD19pCD27mcell1_22_CGAGGCTG.TATCTC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACAACC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACCGCG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACGTGG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ACTCAC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.AGGATG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ATAGCG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.ATCGAC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CAAGAG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CATGAC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CCTTCG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CGGTAG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CTATTG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.CTCAGC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GACACG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GCATTC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GTGAGG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.GTTGAG",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.TATCTC",
                        "RRBS_NormalBCD19pCD27mcell23_44_GTAGAGGA.TCTCTG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ACAACC",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ACGTGG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ACTCAC",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.AGGATG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ATAGCG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.ATCGAC",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CAAGAG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CATGAC",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CCTTCG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CGGTAG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CTATTG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.CTCAGC",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.GACACG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.GTGAGG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.TAGCGG",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.TATCTC",
                        "RRBS_NormalBCD19pCD27mcell45_66_TAAGGCGA.TCTCTG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACAACC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACCGCG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACGTGG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ACTCAC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.AGGATG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ATAGCG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.ATCGAC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CAAGAG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CATGAC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CCTTCG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CGGTAG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CTATTG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.CTCAGC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GACACG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GCATTC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GGCATC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GTGAGG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.GTTGAG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.TAGCGG",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.TATCTC",
                        "RRBS_NormalBCD19pCD27mcell67_88_CGTACTAG.TCTCTG",
                        "RRBS_NormalBCD19pCD27pcell1_22_TAGGCATG.ACAACC",
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
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.GTTGAG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TAGCGG",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TATCTC",
                        "RRBS_NormalBCD19pcell1_22_TAAGGCGA.TCTCTG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACAACC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACCGCG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACGTGG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ACTCAC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.AGGATG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ATAGCG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.ATCGAC",
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
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TAGCGG",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TATCTC",
                        "RRBS_NormalBCD19pcell23_44_CGTACTAG.TCTCTG",
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
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GCATTC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GCTGCC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GGCATC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GTGAGG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.GTTGAG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TAGCGG",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TATCTC",
                        "RRBS_NormalBCD19pcell67_88_TCCTGAGC.TCTCTG",
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
                        'RRBS_cw154_Tris_protease_GR_CAGAGAGG.TCTCTG',
                        'RRBS_trito_pool_1_TAAGGCGA.ACAACC',
                        'RRBS_trito_pool_1_TAAGGCGA.ACGTGG',
                        'RRBS_trito_pool_1_TAAGGCGA.ACTCAC',
                        'RRBS_trito_pool_1_TAAGGCGA.ATAGCG',
                        'RRBS_trito_pool_1_TAAGGCGA.ATCGAC',
                        'RRBS_trito_pool_1_TAAGGCGA.CAAGAG',
                        'RRBS_trito_pool_1_TAAGGCGA.CATGAC',
                        'RRBS_trito_pool_1_TAAGGCGA.CCTTCG',
                        'RRBS_trito_pool_1_TAAGGCGA.CGGTAG',
                        'RRBS_trito_pool_1_TAAGGCGA.CTATTG',
                        'RRBS_trito_pool_1_TAAGGCGA.GACACG',
                        'RRBS_trito_pool_1_TAAGGCGA.GCATTC',
                        'RRBS_trito_pool_1_TAAGGCGA.GCTGCC',
                        'RRBS_trito_pool_1_TAAGGCGA.GGCATC',
                        'RRBS_trito_pool_1_TAAGGCGA.GTGAGG',
                        'RRBS_trito_pool_1_TAAGGCGA.GTTGAG',
                        'RRBS_trito_pool_1_TAAGGCGA.TAGCGG',
                        'RRBS_trito_pool_1_TAAGGCGA.TATCTC',
                        'RRBS_trito_pool_1_TAAGGCGA.TCTCTG',
                        'RRBS_trito_pool_1_TAAGGCGA.TGACAG',
                        'RRBS_trito_pool_1_TAAGGCGA.TGCTGC',
                        'RRBS_trito_pool_2_CGTACTAG.ACAACC',
                        'RRBS_trito_pool_2_CGTACTAG.ACGTGG',
                        'RRBS_trito_pool_2_CGTACTAG.ACTCAC',
                        'RRBS_trito_pool_2_CGTACTAG.AGGATG',
                        'RRBS_trito_pool_2_CGTACTAG.ATAGCG',
                        'RRBS_trito_pool_2_CGTACTAG.ATCGAC',
                        'RRBS_trito_pool_2_CGTACTAG.CAAGAG',
                        'RRBS_trito_pool_2_CGTACTAG.CATGAC',
                        'RRBS_trito_pool_2_CGTACTAG.CCTTCG',
                        'RRBS_trito_pool_2_CGTACTAG.CGGTAG',
                        'RRBS_trito_pool_2_CGTACTAG.CTATTG',
                        'RRBS_trito_pool_2_CGTACTAG.GACACG',
                        'RRBS_trito_pool_2_CGTACTAG.GCATTC',
                        'RRBS_trito_pool_2_CGTACTAG.GCTGCC',
                        'RRBS_trito_pool_2_CGTACTAG.GGCATC',
                        'RRBS_trito_pool_2_CGTACTAG.GTGAGG',
                        'RRBS_trito_pool_2_CGTACTAG.GTTGAG',
                        'RRBS_trito_pool_2_CGTACTAG.TAGCGG',
                        'RRBS_trito_pool_2_CGTACTAG.TATCTC',
                        'RRBS_trito_pool_2_CGTACTAG.TCTCTG',
                        'RRBS_trito_pool_2_CGTACTAG.TGACAG']
      
print(total_matrix.shape)



total_matrix = total_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))



total_matrix = total_matrix.astype(str).apply(''.join)




tott = pd.Series(total_matrix.index.astype(str).str.cat(total_matrix.astype(str),'    '))


tott.to_csv("total_chrom17.phy", header=None, index=None)

print(tott.shape)

             

