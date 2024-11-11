#include "biocma_cst_config.hpp"
#include "common/execinfo.hpp"
#include "core/case_data.hpp"
#include "core/global_initaliser.hpp"
#include "core/simulation_parameters.hpp"

#include <Kokkos_Core.hpp>
#include <api/api.hpp>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>
constexpr int ID_VERIF = 2025;

#ifdef NO_MPI
namespace MPI_W
{
  bool is_initialized() noexcept
  {
    return false;
  }
} // namespace MPI_W
#else
#  include "mpi_w/message_t.hpp"
#endif

// namespace Api
// {

Handle::Handle(uint32_t n_rank, uint32_t current_rank, uint64_t id,uint32_t thread_per_process) 
{
  _data.exec_info = {
      .n_rank = n_rank, .current_rank = current_rank, .thread_per_process = thread_per_process, .verbose = true, .run_id = id};
  id = ID_VERIF;
  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    if (!MPI_W::is_initialized())
    {
      throw std::invalid_argument("INIT MPI FIRST");
    }
  }

  if (!Kokkos::is_initialized())
  {
    Kokkos::initialize(Kokkos::InitializationSettings() 
                           .set_disable_warnings(false)
                           .set_num_threads(static_cast<int32_t>( _data.exec_info.thread_per_process))
                           .set_map_device_id_by("random"));
  }
}


std::optional<std::unique_ptr<Handle>> Handle::init(uint32_t n_rank, uint32_t current_rank, uint64_t id,uint32_t thread_per_proces)
{
  auto ptr = std::unique_ptr<Handle>(new Handle(n_rank, current_rank, id,thread_per_proces));
  return ptr;
}

std::optional<std::unique_ptr<Handle>> Handle::load_mock(uint32_t n_rank, uint32_t current_rank)
{
  auto handle = init(n_rank, current_rank, 0,1).value(); // FIXME
  auto exec = handle->_data.exec_info;

  Core::UserControlParameters mock_param = Core::UserControlParameters::m_default();
  mock_param.final_time = 1000;
  mock_param.delta_time = 0.1;
  mock_param.results_file_name = "./test_api_load";
  mock_param.number_exported_result = 5;
  mock_param.cma_case_path = "/home/benjamin/Documents/code/cpp/BioCMA-MCST/cma_data/0d_mono/";
  mock_param.serde_file = "/home/benjamin/Documents/code/cpp/BioCMA-MCST/results/allo/allo_serde_0.raw";

  auto opt_case = Core::load(exec, std::move(mock_param));

  if (opt_case)
  {
    handle->_data = std::move(*opt_case);
    handle->loaded = true;
    return handle;
  }

  return std::nullopt;
}

Handle::~Handle()
{
  std::cout << "destructing from c++" << std::endl;
}

bool Handle::exec()
{
  if (loaded || (registered && applied))
  {
    try
    {
      std::cout << "Running " << this->_data.exec_info.current_rank << "..." << std::endl;
      Core::exec(std::forward<Core::CaseData>(this->_data));
      return true;
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      return false;
    }
  }
  else
  {
    return false;
  }
}

void Handle::apply(bool to_load)
{
  if (loaded)
  {
    return;
  }

  if (to_load)
  {
    auto opt_case = Core::load(this->_data.exec_info, std::move(this->params));
    if (opt_case)
    {
      this->_data = std::move(*opt_case);
      this->loaded = true;
    }
    return;
  }

  if (!registered)
  {
    throw std::runtime_error("Register first");
    return;
  }

  Core::GlobalInitialiser gi(_data.exec_info, params);
  auto t = gi.init_transitionner();

  auto __simulation = gi.init_simulation();
  if ((!t.has_value() && !__simulation.has_value()) || !gi.check_init_terminate())
  {
    throw std::runtime_error("Error apply");
  }
  _data.params = gi.get_parameters();
  _data.simulation = std::move(*__simulation);
  _data.transitioner = std::move(*t);
  applied = true;
}

void Handle::register_parameters()
{
  if (loaded)
  {
    return;
  }
  this->params = Core::UserControlParameters::m_default();
  this->params.final_time = 10;
  this->params.initialiser_path = "";
  this->params.model_name = "None";
  this->params.delta_time = 0.1;
  this->params.number_particle = 100;
  this->params.results_file_name = "test_api_new";
  this->params.number_exported_result = 50;
  this->params.cma_case_path = "/home/benjamin/Documents/code/cpp/BioCMA-MCST/cma_data/0d_mono/";
  this->params.serde = false;
  // mock_param.flow_files = {mock_param.user_params.cma_case_path};

  registered = true;
}

void Handle::register_parameters(Core::UserControlParameters &&_params)
{
  params = std::move(_params);
  registered = true;
}

bool Handle::register_result_path(std::string_view path)
{
  // TODO Check path
  this->params.results_file_name = path;
  return true; // TODO
}

bool Handle::register_cma_path(std::string_view path)
{
  // TODO Check path
  this->params.results_file_name = path;
  return true; // TODO
}

// } //namespace Api
