#include <cassert>
#include <common/results.hpp>
#include <stdexcept>
#include <string>

using Res = Result<Success, std::string>;

void
test_default_is_valid()
{
  Res r;
  assert(r.valid());
  assert(!r.invalid());
  assert(static_cast<bool>(r));
  // get() on sucess should return default (default string is empty)
  assert(r.get().empty());
}

void
test_error_is_invalid()
{
  Res r{ std::string("boom") };
  assert(!r.valid());
  assert(r.invalid());
  assert(!static_cast<bool>(r));
  assert(r.get() == "boom");
}

void
test_gets()
{
  Res ok;
  // Valid result is ok
  (void)ok.gets();

  // Erreur should throw excpetion
  Res err{ std::string("boom") };
  bool thrown = false;
  try
  {
    (void)err.gets();
  }
  catch (const std::runtime_error&)
  {
    thrown = true;
  }
  assert(thrown && "gets() must throw when invalid");
}

// Rust-like match pattern
void
test_match()
{
  Res ok;
  const int from_ok = ok.match([](auto) { return 1; }, [](auto) { return -1; });
  assert(from_ok == 1);

  Res err{ std::string("boom") };
  const int from_err
      = err.match([](auto) { return 1; }, [](auto) { return -1; });
  assert(from_err == -1);

  // The error branch receives the payload
  Res err2{ std::string("payload") };
  const auto size = err2.match([](auto) { return std::size_t{ 0 }; },
                               [](auto e) { return e.size(); });
  assert(size == std::string("payload").size());
}

int
main()
{
  test_default_is_valid();
  test_error_is_invalid();
  test_gets();
  test_match();
}
