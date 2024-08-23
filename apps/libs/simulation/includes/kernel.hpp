// #ifndef __SIMULATION_KERNEL_HPP__
// #define __SIMULATION_KERNEL_HPP__

// #include "cma_read/light_2d_view.hpp"

// #include "common/common.hpp"
// #include "mc/events.hpp"
// #include "mc/thread_private_data.hpp"
// #include "models/types.hpp"
// #include "transport.hpp"
// #include <Eigen/Core>
// #include <Kokkos_Core.hpp>
// #include <Kokkos_Printf.hpp>
// #include <Kokkos_Random.hpp>
// #include <cstdio>
// #include <mc/container_state.hpp>
// #include <mc/particles/mcparticles.hpp>

// static const std::vector<size_t> v_tmp_index_leaving_flow = {10};
// static const std::vector<double> v_tmp_leaving_flow = {
//     0.}; //{0.000011758 / 10};

// struct Kernel
// {

//   using value_type = MC::Results;

//   std::span<MC::Particles> to_process;
//   KModel kmodel;
//   Kokkos::Random_XorShift1024_Pool<> random_pool;
//   double d_t{};
//   std::span<const double> diag;
//   CmaRead::L2DView<const double> view_cumulative_probability;
//   Kokkos::View<MC::ContainerState *, Kokkos::LayoutRight> domain;
//   Kokkos::View<const size_t **, Kokkos::LayoutStride> neighbors_view;
//   Kokkos::View<double **, Kokkos::LayoutLeft> contribs;

//   Kernel(auto model, auto random_pool, auto _domain_view)
//       : kmodel(model), random_pool(random_pool),
//         domain(std::move(_domain_view)){};

//   ~Kernel() = default;
//   Kernel() = default;
//   Kernel(const Kernel &) = default;
//   Kernel(Kernel &&) noexcept = default;
//   Kernel &operator=(const Kernel &) = default;
//   Kernel &operator=(Kernel &&) noexcept = default;

//   void update(auto _to_process,
//               double _d_t,
//               auto _diag,
//               auto _view_cumulative_probability,
//               auto _neighbors_view,
//               auto _contribs)
//   {
//     to_process = _to_process;
//     d_t = _d_t;
//     diag = _diag;
//     view_cumulative_probability = _view_cumulative_probability;
//     neighbors_view = _neighbors_view;
//     contribs = _contribs;
//   }

//   /**
//    * @brief Updates the state of a particle by moving it to a new compartment
//    * based on probabilistic transitions.
//    *
//    * This function performs a move operation on a particle within a
//    * compartmentalized system. It determines whether the particle should leave
//    * its current compartment based on a probability check and updates its
//    * position if necessary. The new compartment is determined based on a random
//    * number and transition probabilities.
//    *
//    * @tparam T Type of the particle.
//    * @param random_number Random number used for deciding if the particle leaves
//    * the current compartment.
//    * @param random_number2 Random number used for determining the next
//    * compartment.
//    * @param particle Reference to the particle being moved.
//    * @param d_t Time step for the simulation.
//    * @param diag_transition Array of diagonal transition probabilities for each
//    * compartment.
//    * @param cumulative_probability Cumulative probability view used to determine
//    * the next compartment.
//    * @param domain_view Kokkos view of the current state of each compartment.
//    * @param view_neighbors Kokkos view of the neighbors for each compartment.
//    */
//   template <typename DataType>
//   KOKKOS_FUNCTION void
//   sub_kernel_move(double random_number,
//                   double random_number2,
//                   MC::BaseParticle<DataType> &particle) const;

//   template <typename DataType>
//   KOKKOS_FUNCTION void
//   sub_kernel_exit(double random_number,
//                   MC::BaseParticle<DataType> &particle) const;

//   KOKKOS_FUNCTION
//   void operator()(auto &&i_particle, MC::Results &results) const
//   {

//     auto &particle = to_process[i_particle];
//     auto &status = particle.status;
//     if (status == MC::CellStatus::DEAD)
//     {
//       return;
//     }

//     Kokkos::View<double *> vlocal_contribs("local_contribs",
//                                            contribs.extent(0));
//     Kokkos::deep_copy(vlocal_contribs, 0);

//     Eigen::MatrixXd local_contribs = Eigen::Map<Eigen::MatrixXd>(
//         vlocal_contribs.data(), EIGEN_INDEX(contribs.extent(0)), 1);



//     auto generator = random_pool.get_state();
//     const double random_number_1 = generator.drand(0., 1.);
//     const double random_number_2 = generator.drand(0., 1.);
//     const double random_number_3 = generator.drand(0., 1.);
//     random_pool.free_state(generator);

//     sub_kernel_move(random_number_1, random_number_2, particle);

//     sub_kernel_exit(random_number_3, particle);

//     if (status == MC::CellStatus::OUT)
//     {
//       results.events.incr<MC::EventType::Exit>();
//       __ATOM_DECR__(domain[particle.current_container].n_cells)
//       particle.clearState(MC::CellStatus::DEAD);
//       results.index_in_dead_state.emplace_back(i_particle);
//       return;
//     }

//     kmodel.update_kernel(
//         d_t, particle, domain[particle.current_container].concentrations);

//     if (status == MC::CellStatus::DEAD)
//     {
//       results.events.incr<MC::EventType::Death>();
//       __ATOM_DECR__(domain[particle.current_container].n_cells)
//       particle.clearState(MC::CellStatus::DEAD);
//       return;
//     }

//     if (status == MC::CellStatus::CYTOKINESIS)
//     {
//       results.events.incr<MC::EventType::NewParticle>();
//       status = MC::CellStatus::IDLE;

//       __ATOM_INCR__(domain[particle.current_container].n_cells);
//       auto child = kmodel.division_kernel(particle);
//       kmodel.contribution_kernel(child, local_contribs);

//       results.extra_process.emplace_back(std::move(child));
//     }

//     kmodel.contribution_kernel(particle, local_contribs);

//     for (size_t i_species = 0; i_species < contribs.extent(0); ++i_species)
//     {
//       Kokkos::atomic_add(&contribs(i_species, particle.current_container),
//                          local_contribs(EIGEN_INDEX(i_species), 0));
//     }
//   }

//   KOKKOS_INLINE_FUNCTION
//   void join(MC::Results &dst, const MC::Results &src) const
//   {
//     dst.events.inplace_reduce(src.events);

//     dst.extra_process.insert(dst.extra_process.end(),
//                              src.extra_process.begin(),
//                              src.extra_process.end());
//     dst.index_in_dead_state.insert(dst.index_in_dead_state.end(),
//                                    src.index_in_dead_state.begin(),
//                                    src.index_in_dead_state.end());
//   }

//   // KOKKOS_INLINE_FUNCTION
//   // void init(MC::Results &val) const
//   // {
//   //   MC::init_r(val);
//   // }
// };

// template <typename DataType>
// KOKKOS_FUNCTION void
// Kernel::sub_kernel_move(double random_number,
//                         double random_number2,
//                         MC::BaseParticle<DataType> &particle) const

// {
//   // Get the current compartment index and corresponding container state
//   const size_t i_compartment = particle.current_container;
//   const int rowId = static_cast<int>(i_compartment);
//   auto &current_container =
//       domain(i_compartment); // Reference to the current container state

//   // Retrieve the list of neighboring compartments for the current compartment
//   const auto i_neighbor = std::span<const size_t>(
//       Kokkos::subview(neighbors_view, i_compartment, Kokkos::ALL).data(),
//       neighbors_view.extent(1));

//   // Check if the particle will leave the current compartment
//   if (!probability_leaving(random_number,
//                            current_container.volume_liq,
//                            diag[i_compartment],
//                            d_t))
//   {
//     return; // Particle remains in the current compartment
//   }

//   // Find the next compartment based on the random number and cumulative
//   // probabilities
//   const size_t next = Simulation::find_next_compartment(
//       rowId, random_number2, i_neighbor, view_cumulative_probability);

//   // Atomically update the cell counts in the current and next compartments
//   __ATOM_DECR__(
//       current_container
//           .n_cells); // Decrement cell count in the current compartment
//   __ATOM_INCR__(
//       domain(next).n_cells); // Increment cell count in the next compartment

//   // Update the particle's current compartment
//   particle.current_container = next;
// }

// template <typename DataType>
// KOKKOS_FUNCTION void
// Kernel::sub_kernel_exit(double random_number,
//                         MC::BaseParticle<DataType> &particle) const
// {
//   const auto lambda = [=, this, &particle](auto &&index, auto &&flow)
//   {
//     if (particle.current_container != index ||
//         particle.status != MC::CellStatus::IDLE)
//     {
//       return;
//     }

//     if (probability_leaving(random_number, domain[index].volume_liq, flow, d_t))
//     {
//       particle.status = MC::CellStatus::OUT;
//     }
//   };

//   for (size_t i = 0; i < v_tmp_index_leaving_flow.size(); ++i)
//   {
//     const auto &index = v_tmp_index_leaving_flow[i];
//     const auto &flow = v_tmp_leaving_flow[i];
//     lambda(index, flow);
//   }
// }

// #endif
