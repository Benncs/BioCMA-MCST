#ifndef __KOKKOS_EIGEN_HPP__
#define __KOKKOS_EIGEN_HPP__

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <type_traits>

namespace KokkosEigen
{

  static constexpr int CompileMatrixSizeEigen = -1;

  /**
   * @brief Maps a Kokkos layout type to its Eigen storage order equivalent.
   *
   * @tparam _layout Kokkos layout type (LayoutLeft or LayoutRight).
   * @return Eigen::ColMajor for LayoutLeft, Eigen::RowMajor for LayoutRight.
   */
  template <typename _layout>
  static consteval auto
  get_eigen_layout()
  {
    if constexpr (std::is_same_v<_layout, Kokkos::LayoutLeft>)
    {
      return Eigen::ColMajor;
    }
    else if constexpr (std::is_same_v<_layout, Kokkos::LayoutRight>)
    {
      return Eigen::RowMajor;
    }
    else
    {
      static_assert(std::is_same_v<_layout, void>, "Unsupported layout");
    }
  }

  /**
   * @brief Wraps a Kokkos 2D View as an Eigen matrix map
   *
   * @tparam data_type Scalar element type
   * @tparam layout   Kokkos memory layout; defaults to LayoutLeft
   * @tparam exec     Kokkos execution space; defaults to void
   * @tparam args     Additional template arguments for the underlying View
   */
  template <typename data_type,
            typename layout = Kokkos::LayoutLeft,
            typename exec = void,
            typename args = void>
  class KokkosEigen2D
  {

    //!< ColMajor or RowMajor
    static constexpr auto eigen_layout = get_eigen_layout<layout>();
    //!< underlying dual view type
    using view_type = Kokkos::DualView<data_type**, layout, exec, args>;
    //!< dynamic Eigen matrix type
    using eigen_mat_type = Eigen::Matrix<data_type,
                                         CompileMatrixSizeEigen,
                                         CompileMatrixSizeEigen,
                                         eigen_layout>;

    //!< dynamic Eigen array type
    using eigen_array_type = Eigen::Array<data_type,
                                          CompileMatrixSizeEigen,
                                          CompileMatrixSizeEigen,
                                          eigen_layout>;

  public:
    //!< map over raw data type
    using eigen_map_type = Eigen::Map<eigen_mat_type>;

  private:
    view_type m_view;           //!< Kokkos DualView storage
    eigen_map_type m_eigen_map; //!< Eigen map pointing to host data

    template <typename m_type>
    auto
    get_map()
    {
      return m_type(
          m_view.view_host().data(), m_view.extent(0), m_view.extent(1));
    }

  public:
    using view_layout = layout;
    using float_type = data_type;
    //!< underlying host view
    using host_view_type = view_type::t_host;
    //!< underlying device view
    using device_view_type = view_type::t_dev;
    //!< underlying const device view
    using const_device_view_type = view_type::t_dev_const;

    /// @brief Allocates Object of shape (n x m) on host and device.
    /// @param label Identifier name for the Kokkos View.
    /// @param n     Number of rows.
    /// @param m     Number of columns.
    KokkosEigen2D(std::string label, std::size_t n, std::size_t m)
        : m_view(label, n, m), m_eigen_map(get_map<eigen_map_type>())
    {
    }

    KokkosEigen2D() : m_view(), m_eigen_map(nullptr, 0, 0)
    {
    }

    KokkosEigen2D& operator=(const KokkosEigen2D&) = delete;

    KokkosEigen2D& operator=(KokkosEigen2D&& rhs) noexcept = delete;

    KokkosEigen2D(KokkosEigen2D&&) = delete;
    KokkosEigen2D(const KokkosEigen2D&) = delete;
    ~KokkosEigen2D() = default;

    eigen_map_type&
    eigen()
    {
      return m_eigen_map;
    }

    const eigen_map_type&
    cst_eigen() const
    {
      return m_eigen_map;
    }

    view_type&
    view()
    {
      return m_view;
    }

    host_view_type
    host_view()
    {
      return m_view.view_host();
    }

    const_device_view_type
    device_view_cst() const
    {
      return m_view.view_device();
    }

    device_view_type
    device_view() const
    {
      return m_view.view_device();
    }

    [[nodiscard]] std::size_t
    n_row() const
    {
      return m_view.extent(0);
    }
    [[nodiscard]] std::size_t
    n_col() const
    {
      return m_view.extent(1);
    }

    void
    host_to_device_sync()
    {
      m_view.template modify<typename view_type::host_mirror_space>();
      m_view.template sync<typename view_type::execution_space>();
    }

    void
    device_to_host_sync()
    {
      m_view.template modify<typename view_type::execution_space>();
      m_view.template sync<typename view_type::host_mirror_space>();
    }

    eigen_array_type
    as_array()
    {
      return m_eigen_map.array();
    }

    eigen_array_type
    as_array() const
    {
      return m_eigen_map.array();
    }

    void
    set_host(std::size_t i_r, std::size_t i_c, data_type&& value)
    {
      m_view.view_host()(i_r, i_c) = value;
    }
    data_type&
    get_host(std::size_t i_r, std::size_t i_c)
    {
      return m_view.view_host()(i_r, i_c);
    }

    [[nodiscard]] std::span<const double>
    get_span() const
    {
      const auto size = n_col() * n_row();
      return { m_view.view_host().data(), size };
    }

    std::span<double>
    get_span()
    {
      const auto size = n_col() * n_row();
      return { m_view.view_host().data(), size };
    }
  };

  namespace Alias
  {

    /** @brief Template type for eigen diagonal object

    @tparam ftype template is the scalar type needs to be floating point type
    */
    template <typename ftype>
      requires(std::is_floating_point_v<ftype>)
    using DiagonalType = Eigen::DiagonalMatrix<ftype, CompileMatrixSizeEigen>;

    /** @brief Template type for eigen matrix object

    @tparam Layout is the disired data layout (Kokkos::LayoutLeft or
    LayoutRight)
    @tparam ftype template is the scalar type needs to be floating point type
    */
    template <typename Layout, typename ftype>
    using MatrixType = Eigen::Matrix<ftype,
                                     CompileMatrixSizeEigen,
                                     CompileMatrixSizeEigen,
                                     get_eigen_layout<Layout>()>;

    /** @brief Template type for colmajor eigen matrix object

    @tparam ftype template is the scalar type needs to be floating point type
    */
    template <typename ftype>
      requires(std::is_floating_point_v<ftype>)
    using ColMajorMatrixtype = MatrixType<Kokkos::LayoutLeft, ftype>;

    /** @brief Template type for  eigen sparse matrix object

    @tparam ftype template is the scalar type needs to be floating point type
    */
    template <typename Layout, typename ftype>
      requires(std::is_floating_point_v<ftype>)
    using SparseMatrixType
        = Eigen::SparseMatrix<ftype, get_eigen_layout<Layout>()>;

  } // namespace Alias

} // namespace KokkosEigen

#endif