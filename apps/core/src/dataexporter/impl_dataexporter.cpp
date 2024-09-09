#include <string>
#ifdef USE_HIGHFIVE
#  include "common/simulation_parameters.hpp"

#  include <dataexporter/impl_dataexporter.hpp>
#  include <mc/events.hpp>

#  include <Eigen/Core>
#  include <nl_types.h>
#  include <string_view>
#  include <tuple>

#  include <Eigen/Dense>
#  include <highfive/H5DataSpace.hpp>
#  include <highfive/H5File.hpp>
#  include <highfive/H5PropertyList.hpp>
#  include <highfive/eigen.hpp>
#  include <highfive/highfive.hpp>

constexpr size_t hdf5_max_compression = 9;

DataExportHighFive::DataExportHighFive(const ExecInfo &info,
                                       const SimulationParameters &params,
                                       std::string_view _filename,
                                       std::tuple<size_t, size_t> dim,
                                       size_t niter,
                                       std::span<size_t> distribution,
                                       double weight)
    : DataExporter(info, params, _filename, dim, niter, distribution, weight)
{
  auto _file = HighFive::File(filename, HighFive::File::Truncate);
  write_initial(_file, info, initial_values);
  write_attributes(_file, metadata);
  prepare();
}

void DataExportHighFive::write_initial(HighFive::File &file,
                                       const ExecInfo &info,
                                       export_initial_kv &md)
{

  for (const auto &kv : md)
  {
    const std::string &attributeName = "initial_parameters/" + kv.first;
    const export_initial_t &value = kv.second;
    std::visit(
        [&](const auto &val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
          file.createDataSet<T>(attributeName, val);
        },
        value);
  }

  file.createDataSet("misc/n_node_thread", info.thread_per_process);
  file.createDataSet("misc/n_rank", info.n_rank);
}

void DataExportHighFive::write_attributes(HighFive::File &file,
                                          export_metadata_kv &md)
{

  for (const auto &kv : md)
  {
    const std::string &attributeName = kv.first;
    const export_metadata_t &value = kv.second;
    std::visit(
        [&](const auto &val)
        {
          using T = std::decay_t<decltype(val)>; // Get the actual type T of
                                                 // the variant alternative
          file.createAttribute<T>(attributeName, val);
        },
        value);
  }
}

void DataExportHighFive::prepare()
{
  HighFive::File file(filename, HighFive::File::ReadWrite);

  auto create_dataset = [&file](const std::string &name,
                                const std::vector<size_t> &dims,
                                const std::vector<size_t> &max_dims,
                                const HighFive::DataType &dtype,
                                const std::vector<hsize_t> &chunk_dims,
                                size_t compression = hdf5_max_compression)
  {
    HighFive::DataSpace dataspace(dims, max_dims);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk_dims));
    props.add(HighFive::Shuffle());
    props.add(HighFive::Deflate(compression));
    return file.createDataSet(name, dataspace, dtype, props);
  };

  auto double_type = HighFive::create_datatype<double>();
  auto size_t_type = HighFive::create_datatype<size_t>();
  create_dataset("records/concentration_liquid",
                 {1, n_col, n_row},
                 {n_iter, n_col, n_row},
                 double_type,
                 {1, n_col, n_row});
  create_dataset("records/liquid_volume",
                 {1, n_col},
                 {n_iter, n_col},
                 double_type,
                 {1, n_col});
  create_dataset("records/gas_volume",
                 {1, n_col},
                 {n_iter, n_col},
                 double_type,
                 {1, n_col});
  create_dataset("records/number_particles",
                 {1, n_col},
                 {n_iter, n_col},
                 size_t_type,
                 {1, n_col});
  create_dataset("records/time", {1}, {n_iter}, double_type, {1});

  file.flush();
}



void DataExportHighFive::append(double t,
                                std::span<double> concentration_liquid,
                                const std::vector<size_t> &distribution,
                                std::span<const double> liquid_volume,
                                std::span<const double> volume_gas)
{
  try
  {
    HighFive::File file(filename, HighFive::File::ReadWrite);

    auto write_dataset =
        [&file]<typename T>(const std::string &name,
                            const std::vector<size_t> &new_size,
                            const std::vector<size_t> &select_start,
                            const std::vector<size_t> &select_size,
                            std::span<T> data)
    {
      auto dataset = file.getDataSet(name);
      dataset.resize(new_size);
      dataset.select(select_start, select_size).write_raw(data.data());
    };

    auto write_dataset_vector =
        [&file]<typename T>(const std::string &name,
                            const std::vector<size_t> &new_size,
                            const std::vector<size_t> &select_start,
                            const std::vector<size_t> &select_size,
                            const std::vector<T> &data)
    {
      auto dataset = file.getDataSet(name);
      dataset.resize(new_size);
      dataset.select(select_start, select_size).write_raw(data.data());
    };

    auto write_dataset_real =
        [&file]<typename T>(const std::string &name,
                            const std::vector<size_t> &new_size,
                            const std::vector<size_t> &select_start,
                            const std::vector<size_t> &select_size,
                            const T *data)
    {
      auto dataset = file.getDataSet(name);
      dataset.resize(new_size);
      dataset.select(select_start, select_size).write_raw(data);
    };

    write_dataset("records/concentration_liquid",
                  {counter + 1, n_col, n_row},
                  {counter, 0, 0},
                  {1, n_col, n_row},
                  concentration_liquid);
    write_dataset("records/liquid_volume",
                  {counter + 1, n_col},
                  {counter, 0},
                  {1, n_col},
                  liquid_volume);
    write_dataset("records/gas_volume",
                  {counter + 1, n_col},
                  {counter, 0},
                  {1, n_col},
                  volume_gas);
    write_dataset_vector("records/number_particles",
                         {counter + 1, n_col},
                         {counter, 0},
                         {1, n_col},
                         distribution);
    write_dataset_real("records/time", {counter + 1}, {counter}, {1}, &t);

    counter++;
  }
  catch (...)
  {
  }
}

void DataExportHighFive::write_final_results(ExportData &data,
                                             std::span<size_t> distribution)

{

  const auto dim = data.dim;
  const auto nr = static_cast<size_t>(get<0>(dim));
  const auto nc = static_cast<size_t>(get<1>(dim));

  HighFive::File file(filename, HighFive::File::ReadWrite);

  file.createDataSet("final_results/number_particles", data.number_particles);

  auto cl =
      Eigen::Map<Eigen::MatrixXd>(data.concentration_liquid.data(), nr, nc);

  file.createDataSet("final_results/concentrations/liquid", cl);

  if (!data.concentration_gas.empty())
  {
    auto cg =
        Eigen::Map<Eigen::MatrixXd>(data.concentration_gas.data(), nr, nc);
    file.createDataSet("final_results/concentrations/gas", cg);
  }

  file.createDataSet(
      "final_results/distribution",
      std::vector<double>(distribution.begin(), distribution.end()));
  file.createDataSet("final_results/events/move",
                     data.event->get<MC::EventType::Move>());
  file.createDataSet("final_results/events/total_division",
                     data.event->get<MC::EventType::NewParticle>());
  file.createDataSet("final_results/events/total_death",
                     data.event->get<MC::EventType::Death>());

  file.createDataSet("final_results/events/total_exit",
                     data.event->get<MC::EventType::Exit>());

  file.createDataSet<std::string>("biological_model/description ",
                                  std::string("Model description"));
}

void DataExportHighFive::write_particle_data(
    Kokkos::View<std::string *, HostSpace> names,
    Kokkos::View<double **, HostSpace> particle_values,
    Kokkos::View<double **, HostSpace> spatial_values,
    const std::string &ds_name)
{
  HighFive::File file(filename, HighFive::File::ReadWrite);
  HighFive::DataSetCreateProps ds_props;
  const size_t n_particles = particle_values.extent(1);
  const auto n_compartments = spatial_values.extent(1);
  const auto data_space_particle = HighFive::DataSpace(n_particles);
  HighFive::DataSetCreateProps _ds_props_spatial;
  _ds_props_spatial.add(HighFive::Chunking(1));
  _ds_props_spatial.add(HighFive::Shuffle());
  _ds_props_spatial.add(HighFive::Deflate(hdf5_max_compression));

  for (size_t i_name = 0; i_name < names.size(); ++i_name)
  {
    auto *ptr_particles =
        Kokkos::subview(particle_values, i_name, Kokkos::ALL).data();

    auto *ptr_spatial =
        Kokkos::subview(spatial_values, i_name, Kokkos::ALL).data();
    {
      auto dataset =
          file.createDataSet<double>(ds_name + "spatial/" + names[i_name],
                                     HighFive::DataSpace(n_compartments),
                                     _ds_props_spatial);
      dataset.write_raw(ptr_spatial);
    }

    {
      auto ds = file.createDataSet<double>(
          ds_name + names[i_name], data_space_particle, ds_props);
      ds.write_raw(ptr_particles);
    }
  }
}




#endif