#!/bin/bash 

type=release
name="$1"
executable=./builddir/$type/apps/cli/biocma_mcst_cli_app

cli_args=$(./tools/cli_formater.py $name)

./builddir/$type/apps/cli/biocma_mcst_cli_app $cli_args  

# np=1000
# final_time=10

# # $executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ 

# # $executable -np $np -d $final_time -f ./cma_data/avg_14/ -nex 1000 -dt 1e-2

# # $executable -np $np -d $final_time -f ../compartment-modelling-tool/out/sanofi/ -nex 100 

# $executable -np $np -d $final_time -r 1 -f /mnt/c/Users/casale/Documents/cfd/flowmap_14/ -nex 1000 -mn model_monod

# $executable -np $np -d $final_time -r 1  -f  ./cma_data/bench_2/ 

# $executable -np $np -d $final_time -f ./cma_data/sanofi/ 

# $executable -np $np -d $final_time -f ./cma_data/bench/ 

# $executable -np $np -d $final_time -f /home/benjamin/Documenti/cpp/BIREM_Project/out/sanofi/ #-dt 0.0004231431780566232

# $executable -np $np -dt 0.5 -d $final_time -f /mnt/c/Users/casale/Documents/code/cpp/biomc/cma_data/0d/  -nex 1200



# $executable -np $np -d $final_time -f  ../BIREM_new/out/ 


