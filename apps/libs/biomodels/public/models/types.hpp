#ifndef __MODELS_TYPES_HPP__
#define __MODELS_TYPES_HPP__

#include <common/cmodel_parameters.hpp>
#include <mc/particles/mcparticles.hpp>
#include <span>

#include <Eigen/Core>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>


//trivial && standard layour -> POD 
template <typename T, typename DataType>
concept BioModel = requires(std::is_trivial<T>,
                            std::is_standard_layout<T>,
                            T a,
                            MC::BaseParticle<DataType> &p,
                            double d_t,
                            std::span<const double> concentrations,
                            Eigen::MatrixXd &contributionMatrix) {
  { a.init(p) } -> std::same_as<void>;
  { a.update(d_t, p, concentrations) } -> std::same_as<void>;
  { a.divisionl(p) } -> std::same_as<MC::BaseParticle<DataType>>;
  { a.contribution(p, contributionMatrix) } -> std::same_as<void>;
};

#ifdef USE_PYTHON_MODULE
#  include <functional>

using ModelInit = std::function<void(MC::Particles &)>;

using ModelUpdate =
    std::function<void(double, MC::Particles &m, std::span<const double>)>;

using ModelDivision = std::function<MC::Particles(MC::Particles &)>;

using ModelContribution =
    std::function<void(MC::Particles &, Eigen::MatrixXd &)>;

#else
// using ModelUpdate = void (*)(double, MC::Particles &, std::span<const double>);

// using ModelDivision = MC::Particles (*)(MC::Particles &);

// using ModelInit = void (*)(MC::Particles &);

// using ModelContribution = void (*)(MC::Particles &, Eigen::MatrixXd &);

// #endif

// // #ifdef DEBUG
// using ModelDebug = std::function<void(MC::Particles &)>;
// inline void defaut_dgb(MC::Particles & /*unused*/) {};
// // #endif

using model_properties_t = std::variant<double, int, std::string>;

using model_properties_detail_t =
    std::unordered_map<std::string, model_properties_t>;

// using ModelGetProperties =
//     std::function<model_properties_detail_t(const MC::Particles &)>;

// inline model_properties_detail_t
// defaut_properties(const MC::Particles & /*unused*/)
// {
//   return {{"description", "model description"}};
// };

// struct KModel
// {
//   ModelInit init_kernel;
//   ModelUpdate update_kernel;
//   ModelDivision division_kernel;
//   ModelContribution contribution_kernel;
//   ModelGetProperties get_properties = defaut_properties;
// };

#endif //__MODELS_TYPES_HPP__