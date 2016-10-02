import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")


mcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27mcell1_22*")
print(len(mcell))
trito = glob.glob("binary_position_RRBS_trito_pool_1*")
print(len(trito))

totalfiles = mcell + trito
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
                        'RRBS_trito_pool_1_TAAGGCGA.TGCTGC']

print(total_matrix.shape)



total_matrix = total_matrix.applymap(lambda x: int(x) if pd.notnull(x) else str("?"))



total_matrix = total_matrix.astype(str).apply(''.join)




tott = pd.Series(total_matrix.index.astype(str).str.cat(total_matrix.astype(str),'    '))



tott.to_csv("normal_cll_9.phy", header=None, index=None)

print(tott.shape)

             

