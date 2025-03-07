
absolute_path=/home-local/casale/Documents/thesis/code/BioCMA-MCST/tools/..

#flowmap=/home_pers/casale/Documents/thesis/cfd/bench/
# arg_flow=-r 0 -f $flowmap


flowmap=/home_pers/casale/Documents/code/biomc/cma_data/0d_gas/
# arg_flow=-f $flowmap

args=" -f $flowmap -d 1.0 -np 10000000 -dt 0.01 -er bench2 -mn model_monod -nt 1 -force 1"

executable_gpu=$absolute_path/builddir/gpu/apps/cli/biocma_mcst_cli_app 

executable_cpu=$absolute_path/builddir/host/apps/cli/biocma_mcst_cli_app 


# mpirun --report-bindings  --bind-to core  -np 1 $executable_gpu  $args : \
#        --map-by slot:PE=1 -np 5 $executable_cpu $args 

# mpirun --report-bindings --bind-to core -np 1 $executable_gpu $args : \
#         --map-by ppr:5:core -np 5 $executable_cpu $args

mpirun --report-bindings --bind-to core -np 1 $executable_gpu $args : \
        --map-by slot:PE=6  -np15 $executable_cpu $args

# mpirun -np 2 \
#     --bind-to core \
#     --map-by ppr:2:socket \
#     numactl --cpunodebind=0 --membind=0 $executable_cpu $args : \
#     $executable_gpu  $args