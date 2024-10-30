
#include "cma_read/flow_iterator.hpp"
#include "cma_read/flowinfo.hpp"
#include "cma_read/flowmap.hpp"
#include "cma_read/light_2d_view.hpp"
#include "cma_read/neighbors.hpp"
#include "cma_read/reactorstate.hpp"
#include <cstddef>
#include <mpi.h>
#include <mpi_w/impl_op.hpp>
#include <mpi_w/iteration_payload.hpp>
#include <mpi_w/message_t.hpp>

namespace MPI_W
{

  void bcst_reactor_state(CmaRead::ReactorState &data, size_t flow_size, size_t rank)
  {
    // Broadcast the pointer to ReactorState
    MPI_W::broadcast(data.n_compartments, 0);
    size_t fs = data.n_compartments * data.n_compartments;
    std::vector<double> liquid_flow(fs);
    size_t maxneighbor = 0;
    std::vector<size_t> neighbor_data;
    if (rank == 0)
    {
      maxneighbor = data.liquid_flow.getViewNeighors().getNCol();
      auto vi = data.liquid_flow.getViewNeighors();
      neighbor_data = std::vector<size_t>(vi.data().begin(), vi.data().end());
    }
    MPI_W::broadcast(maxneighbor, 0);
    MPI_W::broadcast(neighbor_data, 0, rank);
    MPI_W::broadcast(liquid_flow, 0, rank);

    if (rank != 0)
    {
      auto liq = CmaRead::FlowMap(std::move(liquid_flow), flow_size);
      auto ne = CmaRead::Neighbors(std::move(neighbor_data), flow_size, maxneighbor);
      data.liquid_flow = CmaRead::FlowInfo(std::move(liq), std::move(ne));
    }

    MPI_W::broadcast(data.liquidVolume, 0, rank);
    MPI_W::broadcast(data.gasVolume, 0, rank);
  }

  void bcst_iterator(std::unique_ptr<CmaRead::FlowIterator> &iterator, size_t rank)
  {
    size_t flow_size = 0;
    if (rank == 0)
    {

      // Serialize the FlowIterator object
      std::vector<CmaRead::ReactorState> &data = iterator->getdata();
      flow_size = data[0].n_compartments;
      MPI_W::broadcast(flow_size, 0);
      size_t size_it = data.size();

      // Broadcast the size of the data
      MPI_Bcast(&size_it, sizeof(size_t), MPI_BYTE, 0, MPI_COMM_WORLD);
      // Broadcast each ReactorState pointer
      for (auto &state : data)
      {
        bcst_reactor_state(state, flow_size, rank);
      }
    }
    else
    {
      MPI_W::broadcast(flow_size, 0);
      // Receive the size of the data
      size_t size_it = 0;
      MPI_Bcast(&size_it, sizeof(size_t), MPI_BYTE, 0, MPI_COMM_WORLD);
      // Resize the vector to accommodate the received data
      std::vector<CmaRead::ReactorState> data(size_it);
      // Broadcast each ReactorState pointer
      for (auto &state : data)
      {
        bcst_reactor_state(state, flow_size, rank);
      }
      // iterator = std::make_unique<CmaRead::FlowIterator>(std::move(data)); //FIXME
    }
  }

  class FlowIterator
  {
  public:
    explicit FlowIterator(std::span<std::string> folder, const CmtCommons::cma_exported_paths_t &paths, size_t n_count = 1);

    explicit FlowIterator(std::vector<CmaRead::ReactorState> &&_data);
    double MinLiquidResidenceTime();
    [[nodiscard]] size_t size() const;
    CmaRead::ReactorState &get_unchcked_mut(size_t index);
    [[nodiscard]] const CmaRead::ReactorState &get_unchecked(size_t index) const;

  private:
    std::vector<CmaRead::ReactorState> data;
  };

  IterationPayload::IterationPayload(size_t size_flows, size_t volumes)
  {
    this->liquid_flows.resize(size_flows);
    this->liquid_volumes.resize(volumes);
    this->gas_volumes.resize(volumes);
  }

  void HostIterationPayload::send(size_t rank) const
  {
    int rc1 = MPI_W::send_v<double>(liquid_flows, rank, 0, false);

    int rc2 = MPI_W::send_v<double>(liquid_volumes, rank, 1, false);

    int rc3 = MPI_W::send_v<double>(gas_volumes, rank, 2, false);

    int rc4 = MPI_W::send_v(neigbors.data(), rank, 3, true);

    if (rc1 != MPI_SUCCESS || rc2 != MPI_SUCCESS || rc3 != MPI_SUCCESS || rc4 != MPI_SUCCESS)
    {
      MPI_W::critical_error();
    }
  }

  void IterationPayload::recv(size_t source, MPI_Status *status)
  {

    int rc1 = MPI_W::recv_span<double>(liquid_flows, source, status, 0);
    int rc2 = MPI_W::recv_span<double>(liquid_volumes, source, status, 1);
    int rc3 = MPI_W::recv_span<double>(gas_volumes, source, status, 2);

    auto opt = MPI_W::recv_v<size_t>(source, status, 3);
    if (!opt.has_value())
    {
      MPI_W::critical_error();
    }
    raw_neigbors = opt.value();
    auto n_col = raw_neigbors.size() / gas_volumes.size();
    neigbors = CmaRead::Neighbors::Neighbors_const_view_t(raw_neigbors, gas_volumes.size(), n_col, true);

    if (rc1 != MPI_SUCCESS || rc2 != MPI_SUCCESS || rc3 != MPI_SUCCESS)
    {
      MPI_W::critical_error();
    }
  }

} // namespace MPI_W