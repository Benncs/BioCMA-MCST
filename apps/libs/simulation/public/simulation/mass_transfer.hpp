#ifndef __SIMULATION_HYDRO_MASS_TRANSFER_HPP__
#define __SIMULATION_HYDRO_MASS_TRANSFER_HPP__

namespace MassTransfer
{
  enum class MTRType
  {
    Flowmap,
    FixedKla,
  };

  struct MtrDescriptor
  {
    MTRType type;
  };

} // namespace MassTransfer

#endif