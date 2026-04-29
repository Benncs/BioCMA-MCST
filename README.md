<div align="left">

  [![License: Apache-2.0](https://img.shields.io/badge/License-Apache-blue.svg)](LICENSE)
  [![Version: 1.1.0](https://img.shields.io/badge/Version-1.1.0-red.svg)](LICENSE)
</div>


# BioCMA-MCST: A Biological Simulation Tool Using Monte Carlo and Compartment Modeling Approaches


The **Biological by Compartment Modeling Approach Monte Carlo Simulation Tool** (BioCMA-MCST) aims to provide a efficient simulation tool for chemical engeenireing and bioprocesses for real-world bioreactors.

The main approach is Hybrid Particl-based **Monte Carlo** simulation, which provide framework to solve Population Balance equation applied to biological population. Furthermore, a **Compartment Modelling approach** for spatial resolution provides a good trade-off between CFD of real reactor and integrated/0D/simlplified geometry modelling for details can be found [here](https://github.com/Benncs/rcmtool).

The final goal of this tool is to provide an efficient and user-friendly solution for end-users in the field of biological simulation.

## Project Overview

Simulating biological processes within industrial reactors is computationally demanding: combining the complexity of cellular biochemistry with Computational Fluid Dynamics (CFD) can make costs grow exponentially. Most existing tools address this by either simplifying the fluid dynamics or averaging biological population behavior — trade-offs that risk overlooking critical factors relevant to real industrial settings.
BioCMA-MCST tackles this by decoupling the two challenges. Compartment Modeling allows CFD to be run once and reused across simulations, drastically cutting computational overhead. Meanwhile, the Monte Carlo method captures the natural diversity of biological populations and reactions, producing a more realistic and sensitive model than averaged approaches. Together, these methods deliver detailed, accurate simulations without sacrificing the fidelity needed to understand complex bioprocesses at scale.

## Repository Structure

The repository is organized as follows:
- **apps/**: Contains the main source code for the simulation tool.
- **devutils/**: Includes bootstrap/configuration files.
- **tools/**: Houses various tools that can be used to interface with the code.
- **documentation/**: Contains documentation for both the code and the models used in the simulation.

## Authors

- **CASALE Benjamin**

## Related work

This project is part a PhD work.

If you use this software, please cite it as the following [citation](./CITATION)


### LICENSE

This tool is under [Apache License, Version 2.0](./LICENSE) *(Apache-2.0)*
