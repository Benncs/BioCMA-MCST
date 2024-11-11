#include <api/api.hpp>
#include <api/api_raw.h>
constexpr int ID_VERIF = 2025;

void apply(Handle* handle, int to_load)
{
  if (handle != nullptr)
  {
    handle->apply(to_load != 0);
  }
}

Handle* init_handle_raw(int n_rank, int current_rank, uint64_t id)
{
  auto opt_handle = Handle::init(n_rank, current_rank, id, 1);
  if (opt_handle.has_value())
  {
    return opt_handle->release();
  }
  return nullptr;
}

Handle* load_handle(int n_rank, int current_rank)
{
  auto opt_handle = Handle::load_mock(n_rank, current_rank);
  if (opt_handle.has_value())
  {
    return opt_handle->release();
  }
  return nullptr;
}

void delete_handle(Handle* handle)
{

  if (handle != nullptr)
  {
    std::cout << "Delete..." << std::endl;

    delete handle;
    handle = nullptr;
  }
}

int exec(Handle* handle)
{
  if (handle != nullptr)
  {
    if (handle->get_id() == ID_VERIF)
    {
      return handle->exec() ? 0 : -1;
    }
    return -1;
  }
  return -1;
}

/*
REGISTER
*/

void register_parameters(Handle* handle)
{
  if (handle != nullptr)
  {
    handle->register_parameters();
  }
}

int register_result_path(Handle* handle, const char* c)
{
  if (handle != nullptr)
  {
    return (handle->register_result_path(c)) ? 0 : -1;
  }
  return -1;
}

int register_cma_path(Handle* handle, const char* c)
{
  if (handle != nullptr)
  {
    return (handle->register_cma_path(c)) ? 0 : -1;
  }
  return -1;
}
