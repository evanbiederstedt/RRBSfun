#!/bin/bash

cd /gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/second_test_after_job_failures/standard-RAxML-master/

/gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/second_test_after_job_failures/standard-RAxML-master/raxmlHPC-PTHREADS -T 30 -p 373265 -b 354321 -# 100 -m BINGAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/total_CLL_chrom03.phy -n Tcll_chrom3_multi1
