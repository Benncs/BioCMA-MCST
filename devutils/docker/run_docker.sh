
 cfd_path=../../thesis/cfd-cma/cma_data/
 initaliser_path=./cma_data
 case_path=./tools/cases.xml
 res_path=./results
 docker run -v $initaliser_path:/opt/biomc/cma_data:Z \
 -v $cfd_path:/opt/biomc/ccma_data:Z\
 -v $case_path:/opt/biomc/cases.xml:Z \
 -v $res_path:/opt/biomc/results:Z\
 biocma_omp docker -n 6
