
absolute_path=/home-local/casale/Documents/thesis/code/BioCMA-MCST/tools/..

flowmap=/home_pers/casale/Documents/thesis/cfd/bench/
# arg_flow=-r 0 -f $flowmap


# flowmap=/home_pers/casale/Documents/code/biomc/cma_data/0d_gas/
# arg_flow=-f $flowmap

args="-r 1 -f $flowmap -d 1.0 -np 6000000 -dt 0.001 -er bench2 -mn model_monod -nt 1 -force 1"

executable_gpu=$absolute_path/builddir/gpu/apps/cli/biocma_mcst_cli_app 

executable_cpu=$absolute_path/builddir/host/apps/cli/biocma_mcst_cli_app 




mpirun --report-bindings  --bind-to core:overload-allowed --map-by slot:PE=2  \
         -np 1 $executable_gpu $args  --kokkos-num-threads=1 : \
         -np 2 $executable_cpu $args --kokkos-num-threads=2 : \
         -np 1 $executable_cpu $args --kokkos-num-threads=1


# mpirun --report-bindings  --bind-to core --map-by slot:PE=2  \
#          -np 1 $executable_gpu $args --kokkos_threads=1 : \
#          -np 1 $executable_cpu $args --kokkos_threads=1 : \ 
#          -np 2 $executable_cpu $args --kokkos_threads=2  \ 









