#ifndef __COMMON_HAS_SERIALIZE_HPP__
#define __COMMON_HAS_SERIALIZE_HPP__

#include <type_traits>

namespace
{
  template <typename, typename T> struct __has_serialize
  {
    static_assert(std::integral_constant<T, false>::value,
                  "Second template parameter needs to be of function type.");
  };

  // specialization that does the checking

  template <typename C, typename Ret, typename... Args>
  struct __has_serialize<C, Ret(Args...)>
  {
  private:
    template <typename T>
    static constexpr auto check(T *) ->
        typename std::is_same<
            decltype(std::declval<T>().serde(std::declval<Args>()...)),
            Ret      // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            >::type; // attempt to call it and see if the return type is correct

    template <typename> static constexpr std::false_type check(...);

    using type = decltype(check<C>(0));

  public:
    static constexpr bool value = type::value;
  };
} // namespace
template <class C, typename T> consteval bool has_serialize()
{
  return __has_serialize<C, void(T &)>::value;
}

#endif