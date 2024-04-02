#include <common/common.hpp>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/Sparse>

template <typename T>
using RowMajorDynMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using ColMajorDynMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;



SimulationParameters SimulationParameters::m_default()
{
  SimulationParameters params;
  params.n_particles = 0;
  params.n_species =0;
  params.final_time =0;
  params.flow_files = {"./cma_data/"};
  params.d_t = 0;
  params.n_threads =-1;
  params.verbose = false;
  return params;
}

Eigen::SparseMatrix<double> FlowmapToMat(std::span<double> data, size_t n_row)
{
  if (n_row * n_row != data.size())
  {
    throw std::invalid_argument("FlowMap should be square");
  }
  Eigen::Map<ColMajorDynMatrix<double>> matrix_view(
      data.data(), static_cast<int>(n_row), static_cast<int>(n_row));
  return matrix_view.sparseView();
}
