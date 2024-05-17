#!/bin/bash 

executable=builddir/release/apps/cli/biocma_mcst_cli_app
np=50
final_time=12000

# $executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ 



# $executable -np $np -d $final_time -f ./cma_data/sanofi/ 

# $executable -np $np -d $final_time -f ./cma_data/bioreactor_20m3/ 

$executable -np $np -d $final_time -f /home/benjamin/Documenti/code/cpp/BIREM_new/out/ -dt 1e-2

# $executable -np $np -d $final_time -f  ../BIREM_new/out/ 
