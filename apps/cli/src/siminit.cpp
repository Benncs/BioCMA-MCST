#include <memory>
#include <siminit.hpp>

#include "flow_iterator.hpp"
#include "simulation/models/types.hpp"
#include "simulation/simulation.hpp"
#include <cstdio>
#include <messages/wrap_mpi.hpp>
#include <simulation/models/simple_model.hpp>

#include <mc/mcinit.hpp>
#include <vector>

Simulation::SimulationUnit sim_init(ExecInfo& info,
    SimulationParameters& params,
    std::shared_ptr<FlowIterator>& _flow_handle,
    KModel&& model)
{

    std::vector<double> liq_volume;
    std::vector<double> gas_volume;
    std::vector<std::vector<size_t>> liquid_neighbors;
    size_t n_compartments;
    if (info.current_rank == 0) {
        ReactorState* fstate =nullptr;
        try {
            _flow_handle = std::make_shared<FlowIterator>(params.flow_files);

            size_t n_t = static_cast<size_t>(params.final_time / params.d_t);

            _flow_handle->setRepetition(n_t);

            fstate = _flow_handle->get();
        } catch (const std::exception& e) {
            throw std::runtime_error("Error while reading files");
        }

        n_compartments = fstate->n_compartments;
        liquid_neighbors = _flow_handle->get()->liquid_flow.neigbors;

        liq_volume = fstate->liquidVolume;
        gas_volume = fstate->gasVolume;
    }

    // MPI_W::broadcast(d_t, 0);
    MPI_W::broadcast(liq_volume, 0, info.current_rank);
    MPI_W::broadcast(gas_volume, 0, info.current_rank);

    if (info.current_rank != 0) {
        n_compartments = liq_volume.size();
        liquid_neighbors.resize(n_compartments);
    }

    for (auto&& i : liquid_neighbors) {
        MPI_W::broadcast(i, 0, info.current_rank);
    }

    auto unit = MC::init_unit(info, liq_volume, std::move(liquid_neighbors));

    auto container = MC::init_container(info, params.n_particles);

    auto simulation = Simulation::SimulationUnit(params.n_species,
        std::move(unit),
        std::move(container),
        info,
        info.current_rank == 0);

    simulation.setVolumes(std::move(gas_volume), std::move(liq_volume));

    simulation.post_init(std::move(model));
    return simulation;
}
