
#include <mc/traits.hpp>
#include <species_name_extractor.hpp>
// NOLINTBEGIN
struct ModelWithNames
{
  using uniform_weight = std::true_type; // Using type alias
  using Self = ModelWithNames;
  using FloatType = float;
  using Config = std::nullopt_t;

  enum class particle_var : int
  {
    a = 0,
    b,
    c,
    COUNT
  };

  static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::COUNT);
  static constexpr std::size_t n_c = 1;

  static std::array<std::string, n_c>
  species()
  {
    return { "S" };
  }

  static constexpr std::string_view name = "simple";
  using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;
  using SelfContribs = MC::ParticlesContribs<Self::n_c, Self::FloatType>;

  MODEL_CONSTANT FloatType a_i = 5e-6; // m

  KOKKOS_INLINE_FUNCTION static void
  init([[maybe_unused]] const MC::pool_type& a,
       [[maybe_unused]] std::size_t b,
       [[maybe_unused]] const SelfParticle& c)
  {
  }

  KOKKOS_INLINE_FUNCTION static double
  mass([[maybe_unused]] std::size_t idx,
       [[maybe_unused]] const SelfParticle& arr)
  {
    return 0.;
  }

  KOKKOS_INLINE_FUNCTION static MC::Status
  update([[maybe_unused]] const MC::pool_type& random_pool,
         [[maybe_unused]] FloatType d_t,
         [[maybe_unused]] std::size_t idx,

         [[maybe_unused]] const SelfParticle& arr,
         [[maybe_unused]] const SelfContribs& ac,
         [[maybe_unused]] const std::size_t& position,
         [[maybe_unused]] const MC::LocalConcentration& c)
  {

    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION static void
  division([[maybe_unused]] const MC::pool_type& random_pool,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] std::size_t idx2,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const SelfParticle& buffer_arr)
  {
  }
};
static_assert(ModelType<ModelWithNames>, "Check non ModelWithNames");
static_assert(has_species_name<ModelWithNames>, "Check name of ModelWithNames");

struct ModelWONames
{
  using uniform_weight = std::true_type; // Using type alias
  using Self = ModelWithNames;
  using FloatType = float;
  using Config = std::nullopt_t;

  enum class particle_var : int
  {
    a = 0,
    b,
    c,
    COUNT
  };

  static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::COUNT);
  static constexpr std::size_t n_c = 4;

  static constexpr std::string_view name = "simple";
  using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;
  using SelfContribs = MC::ParticlesContribs<Self::n_c, Self::FloatType>;

  MODEL_CONSTANT FloatType a_i = 5e-6; // m

  KOKKOS_INLINE_FUNCTION static void
  init([[maybe_unused]] const MC::pool_type& a,
       [[maybe_unused]] std::size_t b,
       [[maybe_unused]] const SelfParticle& c)
  {
  }

  KOKKOS_INLINE_FUNCTION static double
  mass([[maybe_unused]] std::size_t idx,
       [[maybe_unused]] const SelfParticle& arr)
  {
    return 0.;
  }

  KOKKOS_INLINE_FUNCTION static MC::Status
  update([[maybe_unused]] const MC::pool_type& random_pool,
         [[maybe_unused]] FloatType d_t,
         [[maybe_unused]] std::size_t idx,

         [[maybe_unused]] const SelfParticle& arr,
         [[maybe_unused]] const SelfContribs& ac,
         [[maybe_unused]] const std::size_t& position,
         [[maybe_unused]] const MC::LocalConcentration& c)
  {

    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION static void
  division([[maybe_unused]] const MC::pool_type& random_pool,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] std::size_t idx2,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const SelfParticle& buffer_arr)
  {
  }
};
static_assert(ModelType<ModelWONames>, "Check non ModelWithNames");
static_assert(!has_species_name<ModelWONames>, "Check name of ModelWithNames");

int
main()
{
  assert(ModelWithNames::species().size() == 1);

  auto list_1 = impl::get_species_names_impl<ModelWithNames>();

  assert(list_1[0] == "S");

  auto list_2 = impl::get_species_names_impl<ModelWONames>();

  assert(list_2.size() == ModelWONames::n_c);
  assert(list_2[0] == "0");
  assert(list_2[list_2.size() - 1] == std::to_string(ModelWONames::n_c - 1));
}
