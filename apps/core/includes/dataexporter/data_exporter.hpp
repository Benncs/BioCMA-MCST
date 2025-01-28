#ifndef __CORE_DATA_EXPORTER_HPP__
#define __CORE_DATA_EXPORTER_HPP__

#include <common/execinfo.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

#include "cmt_common/macro_constructor_assignment.hpp"
#include "traits/Kokkos_IterationPatternTrait.hpp"
#include <span> 
namespace Core
{

  static inline std::string
  get_filename(std::string_view filename,
               std::size_t current_rank)
  {
    return std::string(filename) + "_partial_" + std::to_string(current_rank) + ".h5";
  }

  /**
   * @class DataExporter
   * @brief A class responsible for exporting various types of data such as
   * matrices, arrays, and metadata in a structured format.
   *
   * The `DataExporter` class provides functionality to export numerical data,
   * matrices, arrays, and associated metadata. It supports various data types
   * including numeric, string, and vector-based data, and offers compression
   * and chunking options.
   *
   * The class uses several utility types such as `std::variant` and
   * `std::optional` to handle different types of data. It also includes methods
   * for linking data and organizing it into hierarchies within a file system or
   * storage format.
   *
   */
  class DataExporter
  {
  public:
    DELETE_COPY_MOVE_AC(DataExporter)
    DataExporter() = default;
    /**
     * @brief Creates a link to a specified file.
     *
     * This function can create either a hard link or a symbolic link
     * to the specified file. The type of link created depends on the
     * implementation details. If the original file does not exist,
     * an error will be thrown.
     *
     * @param filename The path to the original file to link to.
     * @param link_name The name of the link that will be created.
     * @param groupname Group that is linked.
     */
    void do_link(std::string_view filename,
                 std::string_view link_name,
                 std::string_view groupname);

  protected:
    /**
     * @struct MultiMatrixDescription
     * @brief Describes the properties of a multi-dimensional matrix for export
     * operations.
     *
     * This structure holds the necessary metadata required to define and manage
     * a matrix, including its dimensions, maximum dimensions, chunking options,
     * and data type characteristics. It is used when preparing or working with
     * matrices within the `DataExporter` class.
     */
    struct MultiMatrixDescription
    {
      std::string name;             ///< Dataset name
      std::vector<size_t> dims;     ///< Dataset dimensions
      std::vector<size_t> max_dims; ///< Expected dataset max dimensions
      std::optional<std::vector<unsigned long long>>
          chunk_dims;   ///< Data chunk along each dimension
      bool compression; ///< Matrix data has to be compressed or not
      bool is_integer;  ///< Matrix data is integer type or floating point
    };

    // Using type definitions with aligned comments
    using export_metadata_t = std::variant<uint64_t, int, std::string>; ///< Metadata types for export
    using export_metadata_kv = std::unordered_map<std::string,
                                                  export_metadata_t>;  ///< Key-value pairs for metadata
    using ViewParticleProperties = Kokkos::View<double **,Kokkos::LayoutRight, HostSpace>; ///< View for particle properties
    using simple_export_t = std::variant<size_t,
                                         std::string,
                                         std::vector<size_t>,
                                         double,uint32_t>; ///< Simple export types
    using export_initial_kv = std::unordered_map<std::string,
                                                 simple_export_t>; ///< Initial export key-value pairs

    using matrix_variant_t = std::variant<std::span<const double>,
                                          std::span<const std::size_t>,
                                          double>; ///< Variant for matrix data types

    // Methods

    explicit DataExporter(
        const ExecInfo &info,
        std::string_view _filename,
        std::optional<export_metadata_t> user_description = std::nullopt);
    ~DataExporter();

    void write_matrix(std::string_view name,
                      std::span<const double> values,
                      bool compress = false);

    void write_matrix(std::string_view name,
                      std::span<const double> values,
                      size_t n_row,
                      size_t n_col,
                      bool compress = false);

    void prepare_matrix(MultiMatrixDescription description);

    void append_matrix(std::string_view name, matrix_variant_t data);

    void append_array(std::string_view name,
                      std::span<const double> data,
                      uint64_t last_size = 0);

    void write_properties(std::optional<std::string> specific_dataspace,
                          const export_metadata_kv &values);

    void write_simple(const export_initial_kv &values, std::string_view root);

    void write_simple(std::string specific_dataspace,
                      const simple_export_t &values);

    export_metadata_kv metadata;
    uint64_t export_counter = 0;

  private:
    std::unordered_map<std::string, MultiMatrixDescription> descriptors;
    class impl;
    std::unique_ptr<impl> pimpl;
  };

} // namespace Core
#endif