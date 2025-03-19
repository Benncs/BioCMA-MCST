\page gettingstarted Getting Started

\tableofcontents

# Getting Started

## Dependencies 
This tools is heavily based on 
- [**Kokkos**](https://kokkos.org/) (\cite carter_edwards_kokkos_2014)
- [**CMtool**](https://gitlab.com/codes_tim/compartment-modelling-tool)  
- [**HighFive**](https://github.com/highfive-devs/highfive)  (\cite devresse_highfive_2024)
- [**Pybind11**](https://pybind11.readthedocs.io/en/stable/index.html) (\cite jakob_pybind11_2017)
- [**Cereal**](http://uscilab.github.io/cereal/) (\cite grant_cereal_2017)


Meson build system as well as Meson have to be installed. 
Dependencies can be handle by Meson buildsystem but user can also use system wide installation system specific configuration. 
Pybind11 dependency ,Cereal as well as HighFive are optionals.  

\note 
    Compiler needs to at least support C++20 standard. 
    Tested for gcc+14 and clang++-18 compiler (with CUDA Version: 12.4)
    Needed Python version >=3.10





## Installation 

[See options for flags](#build-options)

By default all targets are compiled :
- CLI shared : CLI executable without MPI dependency(for shared memory only )
- CLI distributed : CLI executable with MPI dependency (for shared memory only )
- Python Module API: handle_module python api (with MPI dependency only)
- C Raw API: shared libary with and without MPI dependency

### CPU only 
- Downlaod the project from github repo

~~~~~~~~~~~~~bash
meson setup builddir $flags$  
~~~~~~~~~~~~~


### Cuda 

Be sure to be in root project 

~~~~~~~~~~~~~bash
CXX=$(pwd)/devutils/kokkos_wrap/wrap_cxx meson setup builddir $optionalflags$ -Duse_cuda=true -Duse_system_kokkos=true  
~~~~~~~~~~~~~

This command is mandatory to detect the correct compiler, this configuration needs kokkos to be installed system wide for specific GPU. 

**wrap_cxx** script is another wrap around wrapper provided by kokkos to be compatible with Meson. It basically gives the correct compiler to script. 

### Build options 

| Option Name           | Type     | Default Value | Description                                                              |
|-----------------------|----------|---------------|--------------------------------------------------------------------------|
| `use_dynamic_module`   | boolean  | `false`       | Determines if dynamic modules should be used. Defaults to `false`.        |
| `use_cuda`             | boolean  | `false`       | Enables CUDA support if set to `true`. Defaults to `false`.               |
| `use_system_kokkos`    | boolean  | `false`       | Decides whether to use the system-installed Kokkos library. Defaults to `false`. |
| `python_install`       | boolean  | `false`       | Installs Python bindings if set to `true`. Defaults to `false`.           |
| `use_cereal`           | boolean  | `true`        | Enables cereal library support for serialization. Defaults to `true`.     |
| `compile_tools`        | boolean  | `false`       | Enables compilation tools if set to `true`. Defaults to `false`.          |
| `ci_execution`         | boolean  | `false`       | Specifies whether CI execution is enabled. Defaults to `false`.           |



## Fast run 

A python script as well as case example are provided to quickly run simulation from CLI
The cases examples are located in 'tools/cases.xml'. 

A simple run can be exectued with : 
~~~~~~~~~~~~~bash
python3 ./tools/runner.py [case_name] -n [number of thread] [-mpi]
~~~~~~~~~~~~~
More information with 

~~~~~~~~~~~~~bash
python3 ./tools/runner.py --help
~~~~~~~~~~~~~