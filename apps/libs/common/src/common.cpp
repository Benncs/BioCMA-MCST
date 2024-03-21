#include <common/common.hpp>
#include <stdexcept>

#include <Eigen/Sparse>
#include <Eigen/Core>


template<typename T>
using RowMajorDynMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename T>
using ColMajorDynMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;


Eigen::SparseMatrix<double > FlowmapToMat(std::span<double>  data,size_t n_row){
  if(n_row*n_row != data.size())
  {
    throw std::invalid_argument("FlowMap should be square");
  }   
  Eigen::Map<ColMajorDynMatrix<double>> matrix_view(data.data(), n_row, n_row); 
  return matrix_view.sparseView();
}

