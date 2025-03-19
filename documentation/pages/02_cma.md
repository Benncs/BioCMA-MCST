\page Documentation

[TOC]

# Compartment Modelling Approach 

Instead of directly solving hydrodynamics equations, this tool utilizes pre-computed flow fields sourced from various origins and aggregated into a unified compartment formalism.

A compartment represents a coarse mesh element that contains a small volume of fluid or gas where conditions are considered homogeneous. Values are integrated over each volume, preserving spatial heterogeneity while simplifying computations. Both liquid and gas phases are treated as Eulerian scalars transported between compartments using flow field data.

For discrete particles, transport is also based on flow fields but incorporates stochastic methods to handle their movement through the domain.

This documentation focuses on the namespace \ref CmaUtils "CmaUtils" and the functor \ref Simulation::KernelInline::MoveFunctor "MoveFunctor".

## CmaUtils

The **CmaUtils** library serves as the interface between **CMtool** and **Kokkos**. While **CMtool** provides lightweight data accessors, it is not optimized for parallel simulation. Many operations require reading flow data efficiently in a parallel context—this is where **CmaUtils** bridges the gap.

The library wraps **CMtool** and provides an abstraction layer over **Kokkos** to handle data reading and processing from flow maps within the simulation. This design ensures that simulation logic remains clean and focused on computations, while **CmaUtils** manages the complexity of data handling.

### Compartment Operations 

Two main operations have to be performed to CMtool Rawdata in order to used them in simulation kernels. 
- Construct a transition matrix based on in/out flow in each compartment 
- Get the cumulative probabilty destination per compartment 

#### Transition Matrix 

One can reformulate the mass balance over a scalar field y governed by compartmental flow velocities into matrix form.
For each compartment ii, the differential mass balance reads:
\f[
    \frac{dV_{i}y_{i}}{dt}=\sum_{j}y_{j}F(j,i)-y_{i}\sum_{j}y_{j}F(i,j)
\f]
Into 
\f[
    \frac{dYV}{dt}=YM
\f]


Implemented in \ref CmaUtils::get_transition_matrix "get_transition_matrix", this function construct the M matrix. 

Algorithmic Notes:

The implemented method employs a straightforward sparse matrix assembly via triplet lists:
- Off-diagonal terms are populated directly from the flow matrix F(i,j).
- Diagonal terms are computed as the negative sum of outbound flows from each compartment, ensuring local mass conservation.
This structure yields a well-posed linear operator driving the Eulerian scalar transpors. For discrete stochastic transport, only the diagonal is used and extract with \ref CmaUtils::get_diag_transition "this function"

#### Cumulative probability 


This step is performed \ref CmaUtils::get_cumulative_probabilities "here". The result is important to decide where the discrete particles moves. 

The aim is to construct of discrete CDF to have to probability to go in a specific compartment knowning that the cell leaves a compartment.
Let L the event to leave current compartment i, and let \f$I_{j}\f$ the event to go the compartment i we have for each neighbors \f$P(I_{j}|L)\f$.

This function construct  the CDF using the fact that \f[P(I_{j}|L)=\frac{F(i,j)}{\sum_{k} F(i,k)}\f]


### FlowMapTransitionner

Additionally, **CmaUtils** implements the \ref CmaUtils::FlowMapTransitionner "FlowMapTransitionner"—an abstraction designed for managing flow map transitions during the simulation. This component provides seamless access to the correct flow map based on the simulation time, abstracting away the details of data formats and file I/O.
The transtionner supports flexible behaviors, such as:
-   Switching to the next flow map after a set time interval.
-   Looping over available flow maps.
-   Implementing caching mechanisms for improved performance.
-   Interpolation between flowmaps.
This flexibility ensures efficient and accurate simulation of field evolution over time.

