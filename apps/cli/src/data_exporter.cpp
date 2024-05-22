#include "common/simulation_parameters.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5PropertyList.hpp"
#include "models/monod.hpp"
#include "simulation/simulation.hpp"
#include <algorithm>
#include <any>
#include <data_exporter.hpp>
#include <mc/events.hpp>

#include <Eigen/Core>
#include <chrono>
#include <nl_types.h>
#include <tuple>

#include <iomanip>

std::string date_time()
{
  std::stringstream ss;
  auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

void DataExporter::write_final_results(
    const Simulation::SimulationUnit &simulation,
    std::span<size_t> distribution)
{
  ExportData data = {
      simulation.mc_unit->container.to_process.size(),
      simulation.getCliqData(),
      simulation.getCgasData(),
      simulation.mc_unit->ts_events.data(),
      simulation.getDim(),
  };

  write_final_results(data, distribution);
  try
  {
    {
      HighFive::File file(filename, HighFive::File::ReadWrite);
      std::vector<double> mass(simulation.mc_unit->container.to_process.size());
      std::transform(simulation.mc_unit->container.to_process.begin(),
                     simulation.mc_unit->container.to_process.end(),
                     mass.begin(),
                     [](auto &&p)
                     {
                       auto &model = std::any_cast<Monod &>(p.data);
                       return model.l;
                     });
      file.createDataSet("final_results/integrated/mass ", mass);
    }
  }
  catch (...)
  {
  }
}

//////////////////////

#ifdef USE_HIGHFIVE
#  include <Eigen/Dense>
#  include <highfive/eigen.hpp>
#  include <highfive/highfive.hpp>
static void write_attributes(HighFive::File &file)
{

  file.createAttribute<int>("file_version", 1);
  file.createAttribute<std::string>("creation_date", date_time());
  file.createAttribute<std::string>("author", std::string("someone"));
  file.createAttribute<std::string>("description",
                                    std::string("Interesting results"));
}

static void write_initial(HighFive::File &file,
                          ExecInfo &info,
                          SimulationParameters &params,
                          std::span<size_t> distribution)
{

  file.createDataSet("initial_parameters/number_particles_0", params.n_particles);
  file.createDataSet("initial_parameters/number_compartment",
                     params.n_compartments);
  file.createDataSet("initial_parameters/final_time", params.final_time);
  file.createDataSet(
      "initial_parameters/particle_distribution",
      std::vector<size_t>(distribution.begin(), distribution.end()));
  file.createDataSet("initial_parameters/delta_time", params.d_t);
  file.createDataSet("initial_parameters/n_map", params.n_different_maps);
  file.createDataSet("misc/n_node_thread", info.thread_per_process);
  file.createDataSet("misc/n_rank", info.n_rank);
}

DataExporter::DataExporter(ExecInfo &info,
                           SimulationParameters &params,
                           std::string &_filename,
                           std::tuple<size_t, size_t> dim,
                           size_t niter,
                           std::span<size_t> distribution)
    : n_row(get<0>(dim)), n_col(get<1>(dim)), n_iter(niter), filename(_filename)

{
  HighFive::File file(filename, HighFive::File::Truncate);
  write_initial(file, info, params, distribution);
  write_attributes(file);
  prepare();
}

void DataExporter::write_final_results(ExportData &data,
                                       std::span<size_t> distribution)
{
  HighFive::File file(filename, HighFive::File::ReadWrite);
  file.createDataSet("final_results/number_particles", data.number_particles);

  auto dim = data.dim;
  auto nr = static_cast<size_t>(get<0>(dim));
  auto nc = static_cast<size_t>(get<1>(dim));

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

  file.createDataSet("final_results/events/total_division",
                     data.event->get<MC::EventType::NewParticle>());
  file.createDataSet("final_results/events/total_death",
                     data.event->get<MC::EventType::Death>());

  file.createDataSet<std::string>("biological_model/description ",
                                  std::string("Model description"));
}
void DataExporter::prepare()
{
  HighFive::File file(filename, HighFive::File::ReadWrite);

  HighFive::DataSpace dataspace = HighFive::DataSpace(
      {
          1,
          n_col,
          n_row,
      },
      {n_iter, n_col, n_row});

  // Use chunking
  HighFive::DataSetCreateProps props;
  props.add(HighFive::Chunking(std::vector<hsize_t>{1, n_col, n_row}));
  props.add(HighFive::Shuffle());
  props.add(HighFive::Deflate(9));

  // Create the dataset
  HighFive::DataSet dataset =
      file.createDataSet("final_results/concentrations/records",
                         dataspace,
                         HighFive::create_datatype<double>(),
                         props);
                         HighFive::DataSet dataset_2 =
      file.createDataSet("final_results/concentrations/records_distribution",
                         dataspace,
                         HighFive::create_datatype<size_t>(),
                         props);
  file.flush();
}

void DataExporter::append(std::span<double> data,const std::vector<size_t>& distribution)
{
  HighFive::File file(filename, HighFive::File::ReadWrite);
  auto dataset = file.getDataSet("final_results/concentrations/records");
  dataset.resize({counter + 1, n_col, n_row});
  dataset.select({counter, 0, 0}, {1, n_col, n_row}).write_raw(data.data());

  dataset = file.getDataSet("final_results/concentrations/records_distribution");
  dataset.resize({counter + 1, n_col, n_row});
  dataset.select({counter, 0, 0}, {1, n_col, n_row}).write_raw(distribution.data());
  counter++;
}

#else
DataExporter::DataExporter(ExecInfo &info,
                           SimulationParameters &params,
                           std::string &_filename,
                           std::tuple<size_t, size_t> dim,
                           size_t niter,
                           std::span<size_t> distribution)
    : filename(_filename.data())

{
  // NOP
}

void DataExporter::append(std::span<double> data,const std::vector<size_t>& distribution)
{
}

void DataExporter::prepare()
{
}

void DataExporter::write_final_results(ExportData &data,
                                       std::span<size_t> distribution)
{
  // NOP
}
#endif