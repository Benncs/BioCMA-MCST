#!/bin/bash 

executable=builddir/release/apps/cli/biocma_mcst_cli_app
np=10000000
final_time=15

$executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ 
# $executable -np $np -d $final_time -f ./cma_data/sanofi/ 

# $executable -np $np -d $final_time -f ./cma_data/bioreactor_20m3/ 

# $executable -np $np -d $final_time -f ./cma_data/ 

# $executable -np $np -d $final_time -f  ../BIREM_new/out/ 
