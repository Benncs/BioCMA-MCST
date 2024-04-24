#!/bin/bash 

executable=builddir/release/apps/cli/BioCMA-MCST_cli_app
np=10000  
final_time=500

# builddir/release/apps/cli/BioCMA-MCST_cli_app -np 10000 -d 1 -dt 1e-3 -f ./cma_data/raw_6612/ 
$executable -np $np -d $final_time -r 1 -f ./cma_data/test2/ 
# builddir/dynmod/apps/cli/BioCMA-MCST_cli_app -np 1000 -d 1.85 -r 1 -f ./cma_data/test2/ 