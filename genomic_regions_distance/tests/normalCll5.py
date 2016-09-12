import glob
import pandas as pd
import numpy as np

import os
os.chdir("/gpfs/commons/home/biederstedte-934/evan_projects/correct_phylo_files")


mcell = glob.glob("binary_position_RRBS_NormalBCD19pCD27mcell*")
print(len(mcell))
trito = glob.glob("binary_position_RRBS_trito_pool*")
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



tott.to_csv("normal_cll_5.phy", header=None, index=None)

print(tott.shape)

             

