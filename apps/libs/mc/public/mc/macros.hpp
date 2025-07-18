#ifndef __MC_MACRO_HPP__
#define __MC_MACRO_HPP__

#define CHECK_MODEL(name)                                                      \
  static_assert(ModelType<name>, #name);                                       \
  static_assert(FloatingPointType<name::FloatType>, " floatl" #name);

#define MODEL_CONSTANT static constexpr

// Utility to get the index from an enum
#define INDEX_FROM_ENUM(e) static_cast<std::size_t>((e))

// Bounds checking macro using static_extent
#define CHECK_BOUND(__index__, __array_name__)                                 \
  static_assert(INDEX_FROM_ENUM(__index__) < __array_name__.static_extent(1),  \
                "Index out of model bound");

// Main macro that uses bounds checking and array access
#define GET_PROPERTY_FROM(__index__, __array_name__, enum_name)                \
  __array_name__(__index__, INDEX_FROM_ENUM(enum_name))

#define GET_PROPERTY(enum_name) GET_PROPERTY_FROM(idx, arr, enum_name)

#define GET_INDEX(size)                                                        \
  std::size_t idx = (team_handle.league_rank() * team_handle.team_size()) +    \
                    team_handle.team_rank();                                   \
  if (idx >= (size))                                                           \
  {                                                                            \
    return;                                                                    \
  }

#endif