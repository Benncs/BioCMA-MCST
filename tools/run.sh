
absolute_path=/home-local/casale/Documents/thesis/code/BioCMA-MCST/tools/..

#flowmap=/home_pers/casale/Documents/thesis/cfd/bench/
# arg_flow=-r 0 -f $flowmap


flowmap=/home_pers/casale/Documents/code/biomc/cma_data/0d_gas/
# arg_flow=-f $flowmap

args=" -f $flowmap -d 1.0 -np 10000000 -dt 0.01 -er bench2 -mn model_monod -nt 1 -force 1"

exe1=$absolute_path/builddir/gpu/apps/cli/biocma_mcst_cli_app 

exe2=$absolute_path/builddir/host/apps/cli/biocma_mcst_cli_app 


mpirun --report-bindings  --bind-to core  -np 1 $exe1  $args : \
       --map-by slot:PE=1 -np 5 $exe2 $args 