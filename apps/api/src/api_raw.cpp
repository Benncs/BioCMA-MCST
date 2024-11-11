#include <api/api.hpp>
#include <api/api_raw.h>
#include <stdexcept>

void register_parameters(Handle *handle)
{
    if (handle!=nullptr)
    {
        handle->register_parameters();
    }
}
void apply(Handle *handle)
{
    if (handle!=nullptr)
    {
        handle->apply();
    }
}

Handle *init_handle_raw(int n_rank,int current_rank)
{
  auto opt_handle = Handle::init(n_rank,current_rank);
  if (opt_handle.has_value())
  {
    return opt_handle->release();
  }
  return nullptr;
}

Handle* load_handle(int n_rank,int current_rank)
{
    auto opt_handle = Handle::load(n_rank,current_rank);
    if (opt_handle.has_value())
    {
      return opt_handle->release();
    }
    return nullptr;
}

void delete_handle(Handle *handle)
{

  if (handle!=nullptr)
  {
    std::cout << "Delete..." << std::endl;

    delete handle;
    handle = nullptr;
  }
}

int exec(Handle *handle)
{
  if (handle != nullptr)
  {
      if(handle->id==2025)
      {
          return handle->exec() ? 0 : -1;
      }
      return -1;
  }
  return -1;
}
