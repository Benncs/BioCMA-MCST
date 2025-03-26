\page Documentation


[TOC]

# Model declaration

A key feature of BioCMAMC is its flexible model definition capability. The framework is designed to work with any model that adheres to the required interface specifications.

As of version v0.7, models must be known at compile time. Dynamic model definition and loading at runtime are not supported. This limitation primarily arises from the challenge of ensuring code safety and compatibility, particularly when targeting GPU architectures. However, this constraint improves code readability, maintainability, and safety—allowing many potential bugs to be detected during compilation.

For a model to be recognized and executed by BioCMAMC, it must implement specific components:

Refer to the example [example_model.cxx](\ref example_model.cxx "Here") for further details.

## Required Model Subkernels

Each model must define the following subkernels, which are functions executed at different stages of the Monte Carlo particle simulation for each MCparticle:

- init: Executed once before the first simulation time step to initialize the model state.
- update: Called during each iteration to update the state of individual particles.
- contribution: Invoked at each iteration to perform scatter-add operations on the liquid source term.
- division: Triggered during particle division events to handle state changes associated with cell division.

These subkernels form the core computational steps of the model and must be explicitly defined for proper integration with the BioCMAMC framework.



\example example_model.cxx Example of minimum model declaration
