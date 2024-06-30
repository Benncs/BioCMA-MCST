#!/bin/bash 

type=release

executable=./builddir/$type/apps/cli/biocma_mcst_cli_app
np=1000000
final_time=5

# $executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ 

# $executable -np $np -d $final_time -f ./cma_data/avg_14/ -nex 1000 -dt 1e-2

# $executable -np $np -d $final_time -f ../compartment-modelling-tool/out/sanofi/ -nex 100 

$executable -np $np -d $final_time -r 1 -f /mnt/c/Users/casale/Documents/cfd/flowmap_14/ -nex 1000 

# $executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ 

# $executable -np $np -d $final_time -f ./cma_data/sanofi/ 

# $executable -np $np -d $final_time -f ./cma_data/bench/ 

# $executable -np $np -d $final_time -f /home/benjamin/Documenti/cpp/BIREM_Project/out/sanofi/ #-dt 0.0004231431780566232

# $executable -np $np -d $final_time -f /home/benjamin/Documenti/code/cpp/biomc/cma_data/0d/ -dt 0.5




# $executable -np $np -d $final_time -f  ../BIREM_new/out/ 
