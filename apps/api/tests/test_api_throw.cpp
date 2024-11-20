#include "common_test.hpp"
#include <api/api.hpp>
#include <new>

//Overload case where problem occurs in malloc
//Usual new operator throw bad_alloc in this case
//Api use new operator nothrow: that is to say return nullptr instead of throwing excpetion  
void* operator new(std::size_t size, const std::nothrow_t& nothrow_value) noexcept
{
  return nullptr;
}

void test_init_throw()
{
  auto handle = Handle::init(n_rank, i_rank, id, nt);
    
  assert(!handle.has_value());
}

int main()
{

  test_init_throw();
}