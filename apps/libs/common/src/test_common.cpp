#include <common/common.hpp>

#include <cstring>
#include <iostream>

// #include <Eigen/Sparse>
// #include <Eigen/Core>

#define lin_index(i_r, i_c, n_x) (i_r) * (n_x) + (i_c)




// template<typename T>
// using RowMajorDynMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// template<typename T>
// using ColMajorDynMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;


// auto foo(std::span<double> rd)
// {

//     Eigen::Map<ColMajorDynMatrix<double>> matrix_view(rd.data(), 3, 3); 
//     return matrix_view.sparseView();
// } 

void test_common()
{
  std::cout << "from common lib" << std::endl;


  // typedef Eigen::MatrixXd MatrixType;
  // typedef Eigen::Map<MatrixType> MapType;
  // typedef Eigen::Map<const MatrixType> MapTypeConst; 
  // std::vector<double> rd = {1.,0.,3.
  //                           ,4.,5.,6.,
  //                           7.,8.,9.};

  // Eigen::SparseMatrix<double> sm = foo(rd);

  // Eigen::VectorXd x(3);
  // x << 1.0, 2.0, 3.0;
  // std::cout<<sm<<std::endl;

  // auto r = sm.row(1);
  // std::cout<<r<<std::endl;
  // std::cout<<r.coeff(1)<<std::endl;
  // r.coeffRef(1) = -9.0;
  // std::cout<<r<<std::endl;

  // std::cout<<sm<<std::endl;


  // auto res = sm*x;
  // std::cout<<res<<std::endl;

  // auto t = sm.data();
  // auto s = t.size();
  // auto p = t.valuePtr();
  // for(size_t i =0;i<s;++i)
  // {
  //   std::cout<<p[i]<<std::endl;
  // }




}