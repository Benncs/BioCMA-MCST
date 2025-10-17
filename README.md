<div align="left">

  [![License: Apache-2.0](https://img.shields.io/badge/License-Apache-blue.svg)](LICENSE)
  [![Version: 0.9.5](https://img.shields.io/badge/Version-0.9.5-red.svg)](LICENSE)
</div>


# BioCMA-MCST: A Biological Simulation Tool Using Monte Carlo and Compartment Modeling Approaches

The **Biological by Compartment Modeling Approach Monte Carlo Simulation Tool** (BioCMA-MCST) is designed to deliver accurate biological simulations for real-world bioreactor scenarios. This tool provides both time and spatial resolution for complex biological processes.

The main approaches to provides these features are **Monte Carlo** simulation, which provide a simple multi-agent framework for complex phenomena like biological strains behaviour. Furthermore, a **Compartment Modelling approach** for spatial resolution provides a good trade-off between CFD of real reactor and integrated/0D/simlplified geometry modelling.

The final goal of this tool is to provide an efficient and user-friendly solution for end-users in the field of biological simulation.

## Project Overview

Biological simulations face a fundamental challenge: biochemical processes at the cellular level are incredibly complex, numerous, and highly variable. Simulating these processes within a real, complex reactor environment compounds the difficulty, as Computational Fluid Dynamics (CFD) — the standard method for such tasks — is notoriously time-consuming. When attempting to model both biological phenomena and reactor dynamics simultaneously, the computational cost can grow exponentially, making it impractical for large-scale simulations.

To address these challenges, most models and tools for biological simulation employ two primary simplifications. The first simplification involves neglecting some hydrodynamic scales and phenomena, avoiding CFD to manage the biological processes more efficiently. The second simplification simplifies the biological models themselves, focusing only on select phenomena or averaging the behavior of entire populations. While these approaches reduce computational demands, they often result in significant trade-offs, potentially filtering out or overlooking critical factors necessary for understanding specific behaviors in an industrial setting.

The aim of BioCMA-MCST is to maintain the accuracy of Computational Fluid Dynamics (CFD) while incorporating detailed biological models. By using Compartment Modeling, users can simulate CFD once and then reuse the results in subsequent simulations, significantly reducing the computational burden [more information](https://compartment-modelling-tool-codes-tim-1414a41277458b7f47f5759968.gitlab.io/).

Additionally, the Monte Carlo method, while computationally intensive, captures the high diversity in biological populations and reactions, resulting in a more sensitive and realistic model compared to those that rely on averaged biological data. This combination of approaches enables BioCMA-MCST to provide detailed and accurate simulations, offering valuable insights into complex biological processes in industrial settings.

## Repository Structure

The repository is organized as follows:
- **apps/**: Contains the main source code for the simulation tool.
- **devutils/**: Includes bootstrap/configuration files.
- **tools/**: Houses various tools that can be used to interface with the code.
- **documentation/**: Contains documentation for both the code and the models used in the simulation.

## Authors

- **CASALE Benjamin**

## Citation

If you use this software, please cite it as the following [citation] file (./CITATION)


### LICENSE

This tool is under [Apache License, Version 2.0](./LICENSE) *(Apache-2.0)*
