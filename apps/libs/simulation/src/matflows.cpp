#include <simulation/matflows.hpp>

namespace Simulation
{
  MatFlow::MatFlow(const FlowMatrixType &_flows, const FlowMatrixType &_tm)
      : flows(_flows), transition_matrix(_tm)
  {
  }

 

  MatFlow &MatFlow::operator=(MatFlow &&rhs) noexcept
  {
    if (this != &rhs)
    {
      flows = std::move(rhs.flows);
      transition_matrix = std::move(rhs.transition_matrix);
    }
    return *this;
  }

  MatFlow::MatFlow(MatFlow &&other) noexcept
      : flows(std::move(other.flows)),
        transition_matrix(std::move(other.transition_matrix))
  {
  }
}; // namespace Simulation
