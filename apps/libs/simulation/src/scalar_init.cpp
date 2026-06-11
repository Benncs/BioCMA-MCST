#include <Kokkos_CheckUsage.hpp>
#include <optional>
#include <scalar_init.hpp>
#include <scalar_simulation.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation_exception.hpp>

namespace
{
  template <typename T>
  std::vector<T>
  layout_right_to_left(std::span<const T> input, size_t rows, size_t cols)
  {
    if (input.size() != rows * cols)
    {
      throw Simulation::BaseSimulationException(
          Simulation::ErrorCodes::MismatchSize);
    }
    std::vector<T> output(input.size());
    for (size_t i = 0; i < rows; ++i)
    {
      for (size_t j = 0; j < cols; ++j)
      {
        output[j * rows + i] = input[i * cols + j];
      }
    }
    return output;
  }

  template <typename T>
  std::vector<T>
  get_correct_layout_buffer(bool change_layout,
                            std::vector<T>&& input,
                            size_t rows,
                            size_t cols)
  {
    std::vector<T> layout_buffer;
    if (change_layout)
    {
      layout_buffer = layout_right_to_left<T>(input, rows, cols);
    }
    else
    {
      layout_buffer = std::move(input);
    }
    return layout_buffer;
  }

} // namespace

namespace Simulation::impl
{
  void
  post_init_concentration_functor(
      bool is_two_phase_flow,
      ScalarInitializer&& scalar_init,
      const std::shared_ptr<ScalarSimulation>& liquid_scalar,
      const std::shared_ptr<ScalarSimulation>& gas_scalar)
  {

    auto cliq = liquid_scalar->get_concentration();
    std::optional<decltype(cliq)> cgas = std::nullopt; // FIXME

    if (is_two_phase_flow)
    {
      assert(gas_scalar != nullptr);
      cgas = gas_scalar->get_concentration();
      assert(cgas->size() != 0);
    }

    const auto& fv = scalar_init.liquid_f_init.value();

    Kokkos::parallel_for(
        "clear_negs",
        Kokkos::MDRangePolicy<
            decltype(cliq)::execution_space,
            Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
            { 0, 0 }, { cliq.extent(0), cliq.extent(1) }),
        KOKKOS_LAMBDA(int i, int j) {
          cliq(i, j) = fv(i, j);

          if (is_two_phase_flow)
          {
            (*cgas)(i, j) = scalar_init.gas_f_init.value()(i, j);
          }
        });
  }

  void
  post_init_concentration_file(
      bool is_two_phase_flow,
      ScalarInitializer&& scalar_init,
      const std::shared_ptr<ScalarSimulation>& liquid_scalar,
      const std::shared_ptr<ScalarSimulation>& gas_scalar,
      bool change_layout)
  {
    if (!scalar_init.liquid_buffer.has_value())
    {
      throw BaseSimulationException(ErrorCodes::BadInitialiser);
    }
    const std::size_t cols
        = scalar_init.liquid_buffer->size() / scalar_init.n_species;

    // std::vector<double> layout_buffer;
    // if (change_layout)
    // {
    //   layout_buffer = layout_right_to_left<double>(
    //       *scalar_init.liquid_buffer, scalar_init.n_species, cols);
    // }
    // else
    // {
    //   layout_buffer = *scalar_init.liquid_buffer;
    // }
    const auto buffer = get_correct_layout_buffer<double>(
        change_layout,
        std::move(*scalar_init.liquid_buffer),
        scalar_init.n_species,
        cols);

    // Keep this even though get_correct_layout_buffer already may throw
    if (!liquid_scalar->deep_copy_concentration(buffer))
    {
      throw BaseSimulationException(ErrorCodes::MismatchSize);
    }

    if (is_two_phase_flow)
    {
      if (!scalar_init.gas_buffer.has_value())
      {
        throw BaseSimulationException(ErrorCodes::BadInitialiser);
      }
      // if (change_layout)
      // {
      //   layout_buffer = layout_right_to_left<double>(
      //       *scalar_init.gas_buffer, scalar_init.n_species, cols);
      // }
      // else
      // {
      //   layout_buffer = *scalar_init.gas_buffer;
      // }
      const auto buffer = get_correct_layout_buffer<double>(
          change_layout,
          std::move(*scalar_init.gas_buffer),
          scalar_init.n_species,
          cols);

      if (!gas_scalar->deep_copy_concentration(buffer))
      {
        throw BaseSimulationException(ErrorCodes::MismatchSize);
      }
    }
  }
} // namespace Simulation::impl
