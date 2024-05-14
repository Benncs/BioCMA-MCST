#include "cma_read/flowmap.hpp"
#include <transport.hpp>

// void test_valid_input()
// {
//   std::vector<double> data = {1.0,
//                               0.0,
//                               0.0,
//                               0.0,
//                               2.0,
//                               0.0,
//                               0.0,
//                               0.0,
//                               3.0}; // Example square matrix data
//   int n_row = 3; // Square root of the size of the data vector
//   try
//   {
//     auto view = FlowMap::FlowMap_const_view_t(data,n_row);
//     auto sparse_matrix = Simulation::flowmap_to_matrix(data, n_row);
//     // Check if the resulting matrix is correct
//     assert(sparse_matrix.rows() == n_row);
//     assert(sparse_matrix.cols() == n_row);
//     for (int i = 0; i < n_row; ++i)
//     {
//       for (int j = 0; j < n_row; ++j)
//       {
//         if (i == j)
//         {
//           assert(sparse_matrix.coeff(i, j) == data[i * n_row + j]);
//         }
//         else
//         {
//           assert(sparse_matrix.coeff(i, j) == 0.0);
//         }
//       }
//     }
//   }
//   catch (const std::invalid_argument &e)
//   {
//     assert(false);
//   }
// }

// // Test function to check if FlowmapToMat throws an exception for invalid input
// void test_invalid_input()
// {
//   std::vector<double> data = {1.0,
//                               0.0,
//                               0.0, // Non-square matrix data
//                               0.0,
//                               2.0,
//                               0.0};
//   int n_row = 2;
//   try
//   {
//     Simulation::flowmap_to_matrix(data, n_row);
//     assert(false);
//   }
//   catch (const std::invalid_argument &e)
//   {
//     static_assert(true);
//   }
// }

// int main()
// {
//   test_valid_input();
//   test_invalid_input();
//   return 0;
// }