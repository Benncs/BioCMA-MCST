#ifndef __ITERATION_PAYLOAD_HPP__
#define __ITERATION_PAYLOAD_HPP__

#include <cma_read/flow_iterator.hpp>
#include <mpi.h>
#include <span>
#include <vector>
#include <cma_read/neighbors.hpp>

namespace MPI_W
{

  // void bcst_iterator(std::unique_ptr<CmaRead::FlowIterator>& iterator,size_t rank);

  class IterationPayload
  {
  public:
    std::vector<double> liquid_flows;
    std::vector<double> liquid_volumes;
    std::vector<double> gas_volumes;
    std::vector<size_t> raw_neigbors;
    CmaRead::Neighbors::Neighbors_const_view_t neigbors;
    
    explicit IterationPayload(size_t size_flows, size_t volumes);
    // explicit IterationPayload()=delete;
    void recv(size_t source, MPI_Status *status);
  };

  class HostIterationPayload
  {
  public:
    std::span<const double> liquid_flows;
    std::span<const double> liquid_volumes;
    std::span<const double> gas_volumes;
    CmaRead::Neighbors::Neighbors_const_view_t neigbors;
    

    void send(size_t rank) const;
  };

} // namespace MPI_W

#endif //__ITERATION_PAYLOAD_HPP__
