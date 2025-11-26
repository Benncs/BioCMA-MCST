\page cma Documentation



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

This operations has been already studied in \cite morchain_dynamic_2024 
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

## Move kernel


This operations has been already studied in \cite morchain_dynamic_2024 and implementation is \ref Simulation::KernelInline::MoveInfo "here".

This kernel performs the following steps: 
- Compute the probability of leaving the current compartment based on the liquid volume, diagonal transition rate, and time step.
- If the particle moves, draw a random number to select the next compartment using the cumulative transition probabilities.
- Find the neighbor using cumulative probability using a random number.
- Move the particle to the selected compartment; otherwise, it stays in the current one.

The probability of leaving the current compartment follows a discrete Poisson law, where the parameter is the compartment's mean residence time.
\f[
    p_{out} = 1 − \exp(-\frac{\Delta t}{\tau})
\f]
Inverse sampling is used to draw from the discrete Poisson law, allowing a no-branching algorithm. 
The probability of leaving the compartment is \ref Simulation::KernelInline::probability_leaving "implemented" using the regular logarithmic function. 
For small simulation time steps compared to the mean residence time, approximations may also be applied.

The leaving condition is given by:
\f[
\left\{
\begin{array}{l}
u \sim \mathcal{U}(0, 1) \\
\Delta t \cdot F > -\ln(u) \cdot V
\end{array}
\right.
\f]

where:

- \f$u\f$ is a uniform random number in \f$(0, 1)\f$,
- \f$\Delta t\f$ is the simulation time step,
- \f$F\f$ is the leaving flow rate,
- \f$V\f$ is the liquid volume of the compartment.

This Poisson law formalism, based on the cumulative distribution function (CDF) of the leaving probability, is useful because it doesn't depend on the spatial structure or the number of neighboring compartments. It works the same whether a compartment has one neighbor or many. Leaving the reactor is, therefore, treated like moving to any other compartment.

### Leaving reactor  
At each time step, a particle can exit the system through a leaving flow.
The total flow out of a compartment is defined as:\f$F=F_{circulation}+F_{out}\f$. 
Where \f$F=F_{circulation}\f$ comes from transition matrix and flowmap and \f$F=F_{out}\f$ 

The probability of leaving the reactor is computed using the same inverse sampling method from the exponential (Poisson) law, replacing the internal flow with the specific flow leaving the reactor.

