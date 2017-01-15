#!/bin/bash

cd /gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/second_test_after_job_failures/standard-RAxML-master/

/gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/second_test_after_job_failures/standard-RAxML-master/raxmlHPC-PTHREADS -T 30 -p 65732 -b 251433 -# 100 -m BINGAMMA -s /gpfs/commons/home/biederstedte-934/evan_projects/CLL_tests/total_CLL_all.phy -n Tcll18
