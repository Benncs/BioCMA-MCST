// #include "simulation/pc_hydro.hpp"
// #include <transport.hpp>

// #include <Eigen/Dense>
// #include <cassert>

// void test_get_transition_matrix()
// {
//   // Test case 1: 2x2 matrix
//   FlowMatrixType flows_2x2(2, 2);
//   flows_2x2 << 1, 2, /* row */ 3, 4;

//   FlowMatrixType expected_2x2(2, 2);
//   expected_2x2 << -2, 2, 3, -3;

//   assert(Simulation::get_transition_matrix(flows_2x2).isApprox(expected_2x2));

//   // Test case 2: 3x3 matrix
//   FlowMatrixType flows_3x3(3, 3);
//   flows_3x3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

//   FlowMatrixType expected_3x3(3, 3);
//   expected_3x3 << -5, 2, 3, 4, -10, 6, 7, 8, -15;

//   assert(Simulation::get_transition_matrix(flows_3x3).isApprox(expected_3x3));
// }

// #include <iostream>
// void test_get_transition_matrix_integrety()
// {

//   Eigen::MatrixXd flows_2x2 =
//       Eigen::MatrixXd::Random(2, 2); // Generate random 2x2 matrix
//   Eigen::MatrixXd flows_n_n =
//       Eigen::MatrixXd::Random(5, 5); // Generate random 5x5 matrix

//   // Scale the random numbers to be between 0 and 10
//   flows_2x2 = (flows_2x2.array() + 1.0) * 5.0; // Scale to be between 0 and
//   10 flows_n_n = (flows_n_n.array() + 1.0) * 5.0; // Scale to be between 0
//   and 10 auto check = [](auto &f)
//   {
//     std::cout << f << std::endl;
//     auto transition = Simulation::get_transition_matrix(f);
//     std::cout << transition << std::endl;
//     // Check if diagonal elements match the definition
//     for (int i = 0; i < f.rows(); ++i)
//     {
//       double sum_i_j = 0;
//       for (int j = 0; j < f.cols(); ++j)
//       {
//         if (i != j)
//         {
//           sum_i_j += f(i, j);
//         }
//       }
//       sum_i_j = -sum_i_j;

//       std::cerr << transition(i, i) << " " << sum_i_j << std::endl;
//       assert(std::abs(transition(i, i) - sum_i_j) < 1e-5);
//     }
//   };

//   check(flows_2x2);
//   check(flows_n_n);
// }

// int main()
// {
//   test_get_transition_matrix();
//   test_get_transition_matrix_integrety();
// }