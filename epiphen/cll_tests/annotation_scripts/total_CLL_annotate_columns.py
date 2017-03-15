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
    df_list.append(df)


print(len(df_list))

total_matrix = pd.concat([df.set_index("position") for df in df_list], axis=1).reset_index().astype(object)

indexed_matrix = total_matrix ## keep a copy for index of genomic coordinates

total_matrix = total_matrix.drop("index", axis=1)

drop_columns = total_matrix  ## keep copy in order to create 0/1/? matrix such that each character is a column

len(drop_columns.columns)

len(total_matrix.columns)


cell_samples  = ['RRBS_cw154_CutSmart_proteinase_K_TAGGCATG.ACAACC',
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

total_matrix.columns = cell_samples

print(total_matrix.shape)

## >>> print(total_matrix.shape)
## (6336559, 104)

drop_columns = drop_columns.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))

drop_columns = drop_columns.astype(str).apply(''.join)

drop_columns  = drop_columns.reset_index()

total_matrix = total_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))

total_matrix = total_matrix.astype(str).apply(''.join)


os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests")

tott = pd.Series(total_matrix.index.astype(str).str.cat(total_matrix.astype(str),'    '))

tott_drop_columns = pd.Series(drop_columns.index.astype(str).str.cat(total_matrix.astype(str),'    '))  ## [104 rows x 6336566 columns]

print(tott.shape)
print(tott_drop_columns.shape)


df_tott_column_position = tott_drop_columns.apply(lambda x: pd.Series(list(x)))  ## [104 rows x 6336566 columns]
## extra NaN's here



df_tott_column_position_T = df_tott_column_position.T   ## create transpose, and shift on columns [I don't think there is a pandas-efficient way to shift row elements left/right systematically]

for i in range(10):  ## 0 to 9
    df_tott_column_position_T[i]= df_tott_column_position_T[i].shift(2)

for i in range(90):
    j = i + 10      ## 10 to 99
    df_tott_column_position_T[j]= df_tott_column_position_T[j].shift(1)


df_tott_column_position = df_tott_column_position_T.T


df_tott_column_position.drop( df_tott_column_position.columns[[i for i in range(7)]], axis=1, inplace=True) ## drop first 6 columns

### rename columns
indexed_matrixT = indexed_matrix.T
df_tott_column_position.columns = indexed_matrixT.ix[0]


integers_to_sort = df_tott_column_position.columns.to_series().str.extract("([a-z-A-Z]+)(\d*)_(\d+)", expand=True)  # use str.extract to get integers to sort

integers_to_sort[1] = integers_to_sort[1].str.zfill(2)
integers_to_sort[2] = integers_to_sort[2].str.zfill(10)
integers_to_sort["new_coordinates"] = integers_to_sort.apply(lambda x: "{}{}_{}".format(x[0],x[1],x[2]), axis=1)

df_tott_column_position.columns =  integers_to_sort["new_coordinates"]

df_tott_column_position.columns.name = None

df_tott_column_position = df_tott_column_position.sort_index(axis=1)

df_tott_column_position.insert(0, "cell_sample", value = cell_samples)

df_tott_column_position.to_csv("total_CLL_genomic_coordinates_sorted.csv", header=None, index=None)





