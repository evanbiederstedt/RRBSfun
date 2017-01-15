#!/bin/bash

cd /gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/second_test_after_job_failures/standard-RAxML-master/

/gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/second_test_after_job_failures/standard-RAxML-master/raxmlHPC-PTHREADS -T 30 -p 673265 -b 654321 -# 100 -m BINGAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/total_CLL_chrom06.phy -n Tcll_chrom6_multi1
