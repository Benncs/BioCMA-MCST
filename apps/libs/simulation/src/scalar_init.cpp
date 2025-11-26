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
} // namespace

namespace Simulation::impl
{
  void post_init_concentration_functor(
      bool is_two_phase_flow,
      const ScalarInitializer& scalar_init,
      const std::shared_ptr<ScalarSimulation>& liquid_scalar,
      const std::shared_ptr<ScalarSimulation>& gas_scalar)
  {

    auto& cliq = liquid_scalar->get_concentration();
    decltype(&cliq) cgas = nullptr; // FIXME

    if (is_two_phase_flow)
    {
      assert(gas_scalar != nullptr);
      cgas = &gas_scalar->get_concentration();
      assert(cgas->size() != 0);
    }

    const auto& fv = scalar_init.liquid_f_init.value();

    for (decltype(cliq.rows()) i_row = 0; i_row < cliq.rows(); ++i_row)
    {
      for (decltype(cliq.cols()) i_col = 0; i_col < cliq.cols(); ++i_col)
      {
        cliq(i_row, i_col) = fv(i_row, i_col);
        if (is_two_phase_flow)
        {
          (*cgas)(i_row, i_col) = scalar_init.gas_f_init.value()(i_row, i_col);
        }
      }
    }
  }

  void post_init_concentration_file(
      bool is_two_phase_flow,
      const ScalarInitializer& scalar_init,
      const std::shared_ptr<ScalarSimulation>& liquid_scalar,
      const std::shared_ptr<ScalarSimulation>& gas_scalar)
  {
    if (!scalar_init.liquid_buffer.has_value())
    {
      throw SimulationException(ErrorCodes::BadInitialiser);
    }
    const std::size_t cols =
        scalar_init.liquid_buffer->size() / scalar_init.n_species;
    const auto layout_left_buffer = layout_right_to_left<double>(
        *scalar_init.liquid_buffer, scalar_init.n_species, cols);

    if (!liquid_scalar->deep_copy_concentration(layout_left_buffer))
    {

      throw SimulationException(ErrorCodes::MismatchSize);
    }

    if (is_two_phase_flow)
    {
      if (!scalar_init.gas_buffer.has_value())
      {
        throw SimulationException(ErrorCodes::BadInitialiser);
      }
      const auto layout_left_buffer = layout_right_to_left<double>(
          *scalar_init.gas_buffer, scalar_init.n_species, cols);
      if (!gas_scalar->deep_copy_concentration(layout_left_buffer))
      {
        throw SimulationException(ErrorCodes::MismatchSize);
      }
    }
  }
} // namespace Simulation::impl
