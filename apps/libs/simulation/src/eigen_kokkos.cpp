// #include <eigen_kokkos.hpp>

// EigenKokkos::EigenKokkos(std::size_t n_row, std::size_t n_col)
// {
//   eigen_data = Eigen::MatrixXd(n_row, n_col);
//   eigen_data.setZero();

//   host = Kokkos::View<double**, Kokkos::LayoutLeft, HostSpace>(eigen_data.data(), n_row, n_col);
//   compute = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host);
// }

// void EigenKokkos::update_host_to_compute() const
// {
//   Kokkos::deep_copy(compute, host);
// }

// void EigenKokkos::update_compute_to_host() const
// {
//   Kokkos::deep_copy(host, compute);
// }
