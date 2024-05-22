#!/bin/bash 

executable=builddir/apps/cli/biocma_mcst_cli_app
np=5000000
final_time=45

$executable -np $np -d $final_time -r 1  -f  ./cma_data/bench/ 



# $executable -np $np -d $final_time -f ./cma_data/sanofi/ 

# $executable -np $np -d $final_time -f ./cma_data/bioreactor_20m3/ 

# $executable -np $np -d $final_time -f /home/benjamin/Documenti/cpp/BIREM_Project/out/sanofi/ #-dt 0.0004231431780566232

# $executable -np $np -d $final_time -f ./cma_data/0_d/ -dt 1e-2




# $executable -np $np -d $final_time -f  ../BIREM_new/out/ 
