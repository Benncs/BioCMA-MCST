#!/bin/bash 

executable=builddir/release/apps/cli/biocma_mcst_cli_app
np=10000
final_time=100

$executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ -dt 1e-2



# $executable -np $np -d $final_time -f ./cma_data/sanofi/ 

# $executable -np $np -d $final_time -f ./cma_data/bench/ 

# $executable -np $np -d $final_time -f /home/benjamin/Documenti/cpp/BIREM_Project/out/sanofi/ #-dt 0.0004231431780566232

# $executable -np $np -d $final_time -f /home/benjamin/Documenti/code/cpp/biomc/cma_data/0d/ -dt 0.5




# $executable -np $np -d $final_time -f  ../BIREM_new/out/ 
