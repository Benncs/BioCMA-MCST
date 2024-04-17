#ifndef __ITERATION_PAYLOAD_HPP__
#define __ITERATION_PAYLOAD_HPP__

#include <vector>
#include <mpi.h>
#include <span> 

namespace MPI_W {

  class IterationPayload
  {
    public:
      std::vector<double> liquid_flows;
      std::vector<double> liquid_volumes;
      std::vector<double> gas_volumes;

      explicit IterationPayload(size_t size_flows,size_t volumes);
      // explicit IterationPayload()=delete;
      void recv(size_t source,MPI_Status* status);
  };

   class HostIterationPayload
  {
    public:
      std::span<double> liquid_flows;
      std::span<double> liquid_volumes;
      std::span<double> gas_volumes;

      void send(size_t rank)const;
  };


} // namespace MPI_W

#endif //__ITERATION_PAYLOAD_HPP__
