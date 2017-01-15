import glob
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 50) # print all rows


import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")



cw154 = glob.glob("binary_position_RRBS_cw154*")
trito = glob.glob("binary_position_RRBS_trito_pool*")



print(len(cw154))
print(len(trito))

totalfiles = cw154 + trito
print(len(totalfiles))

df_list = []
for file in totalfiles:
    df = pd.read_csv(file)
    df = df.drop("Unnamed: 0", axis=1)
    df["chromosome"] = df["position"].map(lambda x: str(x)[:5])
    df = df[df["chromosome"] == "chr8_"]
    df = df.drop("chromosome", axis=1)
    df_list.append(df)


print(len(df_list))


total_matrix = pd.concat([df.set_index("position") for df in df_list], axis=1).reset_index().astype(object)

total_matrix = total_matrix.drop("index", axis=1)

len(total_matrix.columns)


total_matrix.columns = ["RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACAACC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACCGCG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACGTGG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACTCAC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.AGGATG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ATAGCG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ATCGAC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CAAGAG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CATGAC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CCTTCG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CGGTAG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.CTCAGC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GACACG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GCATTC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GCTGCC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GGCATC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.GTGAGG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.TAGCGG",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.TATCTC",
                        "RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.TCTCTG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACAACC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACCGCG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACGTGG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ACTCAC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.AGGATG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ATAGCG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.ATCGAC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CATGAC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CCTTCG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CGGTAG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CTATTG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.CTCAGC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GACACG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GCATTC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GCTGCC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GGCATC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GTGAGG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.GTTGAG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.TAGCGG",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.TATCTC",
                        "RRBS_cw154_Tris_protease_CTCTCTAC.TCTCTG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACAACC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACCGCG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACGTGG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ACTCAC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.AGGATG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ATAGCG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.ATCGAC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CATGAC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CCTTCG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CGGTAG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CTATTG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.CTCAGC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GACACG",
                        "RBS_cw154_Tris_protease_GR_CAGAGAGG.GCATTC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GCTGCC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GGCATC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GTGAGG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.GTTGAG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.TAGCGG",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.TATCTC",
                        "RRBS_cw154_Tris_protease_GR_CAGAGAGG.TCTCTG",
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


os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests")

tott.to_csv("total_CLL_chrom08.phy", header=None, index=None)

print(tott.shape)
