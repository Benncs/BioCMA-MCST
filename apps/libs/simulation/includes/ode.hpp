// #ifndef __SIMULATION_ODE_HPP__

// #  include <Eigen/Dense>

// namespace ODE
// {

//   struct ButcherTable
//   {

//     Eigen::Matrix3d A;
//     Eigen::Matrix<double, 4, 1> B;

//     ButcherTable();
//   };

//   class ExplicitIntegrator
//   {
//   public:
//     ExplicitIntegrator() = default;

//     Eigen::MatrixXd
//     step(double d_t, const Eigen::MatrixXd &dY, const Eigen::MatrixXd &Y);

//   private:
//     ButcherTable table;
//   };

// } // namespace ODE

// #endif