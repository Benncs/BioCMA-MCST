\page benches Benchmarks

**BioCMAMC** was benchmarked using **four** specific test cases. The test cases are designed to represent typical real-world use scenarios
- **0DB_benchmark**: 0D batch reactor with reaction  
- **3DB_benchmark**: 3D batch reactor with reaction  
- **ODC_benchmark**: 0D chemostat reactor with reaction 
- **0DBE_benchmark**: 0D chemostat reactor with reaction with export results

The first three benchmarks focus exclusively on computational performance without data export to isolate the core calculation cost. These tests were designed with a relatively small number of iterations therefore, I/O overhead is significant under these conditions, it becomes negligible in practical simulations where the number of iterations far exceeds the number of data exports.

To evaluate performance consistently, a specific metric was defined: execution time per 1000 iterations. Given that BioCMAMC uses an explicit integration scheme, the number of iterations is straightforward to determine. Typically, 1000 iterations correspond to 0.1 to 1 second of simulated time, providing a practical and scalable measure of computational efficiency.   

Furthermore, these tests were carried out across 4 different hardware configurations to assess the performance of BioCMAMC under various computing environments:
- CPU Only: Fore shared memory parallelization model (SMP/NUMA)
- GPU Only: Leverages a shared memory parallelization model using GPU resources exclusively (GP-GPU).
- MPI CPU Only: Fore distributed memory system model for CPU-based parallelization.
- Hybrid MPI/CPU/GPU: Supports multi-node configurations, NUMA systems, or environments without dedicated computing devices, combining MPI, CPU, and GPU resources.

These configurations were tested to evaluate how BioCMAMC performs across different memory and parallelization schemes, simulating a variety of real-world computational environments. It is important to note that depending on the specific simulation, some architectures may be more efficient than others. Performing benchmarks helps identify the most optimal configuration for each use case, allowing for better-informed decisions on hardware selection.

## Profiling 

Critical section can be profiled using [Kokkos tools](https://github.com/kokkos/kokkos-tools), this offers a lightweight a simple tool that gives us coarse analysing. 



## Scaling 