#include <api/api.hpp>
#include <new>

// Overload case where problem occurs in malloc
// Usual new operator throw bad_alloc in this case
// Api use new operator nothrow: that is to say return nullptr instead of
// throwing excpetion
void* operator new(std::size_t size,
                   const std::nothrow_t& nothrow_value) noexcept
{
  return nullptr;
}
int main(int argc, char** argv)
{
  auto handle = Api::SimulationInstance::init(argc, argv);
  assert(!handle.has_value());
}