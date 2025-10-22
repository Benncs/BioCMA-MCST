#include "common/logger.hpp"
#include <cli_parser.hpp>
#include <core/simulation_parameters.hpp>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

static CliResults<Core::UserControlParameters>
sanitise_check_cli(Core::UserControlParameters&& params);
static void print_red(std::ostream& os, std::string_view message);
static void throw_bad_arg(std::string_view arg);

static CliResults<Core::UserControlParameters>
parseArg(Core::UserControlParameters& user_controll,
         std::string_view current_param,
         std::string_view current_value);

static CliResults<Core::UserControlParameters> parse_user_param(
    const std::shared_ptr<IO::Logger>& logger, int argc, char** argv)
{
  using return_type = CliResults<Core::UserControlParameters>;
  Core::UserControlParameters control =
      Core::UserControlParameters::m_default();
  if (argc <= 1)
  {
    return return_type("Need at leats one argument");
  }
  try
  {
    int iarg = 1;
    while (iarg < argc - 1)
    {
      std::string_view current_param = argv[iarg];     // NOLINT
      std::string_view current_value = argv[iarg + 1]; // NOLINT
      if (current_param.data() != nullptr && current_param[0] != '\0')
      {
        if (current_param == "h")
        {
          if (logger)
          {
            logger->raw_log(get_help_message());
          }
          exit(0);
        }
        auto opt = parseArg(control, current_param, current_value);
        if (opt)
        {
          control = opt.gets();
        }
        else
        {
          return opt;
        }
      }
      else
      {
        return return_type("ALLO");
      }
      iarg += 2;
    }
    return return_type(std::move(control));
  }
  catch (std::invalid_argument& e)
  {
    std::cerr << e.what() << '\n';
    return return_type(e.what());
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what() << '\n';
    return return_type(e.what());
  }
}

CliResults<Core::UserControlParameters> parse_cli(
    const std::shared_ptr<IO::Logger>& logger, int argc, char** argv) noexcept
{
  auto opt_control = parse_user_param(logger, argc, argv);
  if (!opt_control)
  {
    return CliResults<Core::UserControlParameters>(opt_control.get());
  }
  return sanitise_check_cli(opt_control.gets());
}

static CliResults<Core::UserControlParameters>
parseArg(Core::UserControlParameters& user_control,
         std::string_view _current_param,
         std::string_view current_value)
{
  std::string_view current_param = _current_param.substr(1);

  // Create a map to associate options with corresponding actions (lambdas)
  static const std::unordered_map<std::string_view,
                                  std::function<void(std::string_view)>>
      param_handlers = {
          {"er",
           [&user_control](std::string_view value)
           { user_control.results_file_name = std::string(value); }},
          {"mn",
           [&user_control](std::string_view value)
           { user_control.model_name = std::string(value); }},
          {"nt",
           [&user_control](std::string_view value)
           { user_control.n_thread = std::stoi(std::string(value)); }},
          {"np",
           [&user_control](std::string_view value)
           { user_control.number_particle = std::stol(std::string(value)); }},
          {"nex",
           [&user_control](std::string_view value)
           {
             user_control.number_exported_result =
                 std::stol(std::string(value));
           }},
          {"dt",
           [&user_control](std::string_view value)
           { user_control.delta_time = std::stod(std::string(value)); }},
          {"d",
           [&user_control](std::string_view value)
           { user_control.final_time = std::stod(std::string(value)); }},
          {"r",
           [&user_control](std::string_view)
           { user_control.recursive = true; }},
          {"f",
           [&user_control](std::string_view value)
           { user_control.cma_case_path = std::string(value); }},
          {"force",
           [&user_control](std::string_view)
           { user_control.force_override = true; }},
          {"fi",
           [&user_control](std::string_view value)
           { user_control.initialiser_path = std::string(value); }},
          {"serde",
           [&user_control](std::string_view value)
           {
             user_control.load_serde = true;
             user_control.serde_file = std::string(value);
           }}};

  auto it = param_handlers.find(current_param);
  if (it != param_handlers.end())
  {
    it->second(current_value);
  }
  else
  {
    return CliResults<Core::UserControlParameters>(std::string(current_param));
  }

  return CliResults<Core::UserControlParameters>(std::move(user_control));
}

// static void parseArg(Core::UserControlParameters& user_control,
//                      std::string_view _current_param,
//                      std::string_view current_value)
// {
//   std::string path;
//   // check that begin()+1 is not OOB
//   if (_current_param.size() == 1)
//   {
//     return;
//   }

//   auto current_param = std::string(_current_param.begin() + 1,
//   _current_param.end()); switch (current_param[0])
//   {
//   case 'e':
//   {
//     if (current_param == "er")
//     {
//       user_control.results_file_name = std::string(current_value);
//     }
//     break;
//   }
//   case 'm':
//   {
//     if (current_param == "mn")
//     {
//       user_control.model_name = std::string(current_value);
//     }
//     break;
//   }
//   case 'n':
//   {
//     if (current_param == "nt")
//     {
//       user_control.n_thread = std::stoi(std::string(current_value));
//     }
//     else if (current_param == "np")
//     {
//       user_control.number_particle = std::stol(std::string(current_value));
//     }
//     else if (current_param == "nex")
//     {
//       user_control.number_exported_result =
//       std::stol(std::string(current_value));
//     }
//     break;
//   }
//   case 'd':
//   {
//     if (current_param == "dt")
//     {
//       user_control.delta_time = std::stod(std::string(current_value));
//     }
//     else if (current_param == "d")
//     {
//       user_control.final_time = std::stod(std::string(current_value));
//     }
//     break;
//   }

//   case 'r':
//   {
//     user_control.recursive = true;
//     break;
//   }
//   case 'f':
//   {
//     if (current_param == "f")
//     {
//       user_control.cma_case_path = current_value;
//     }
//     else if (current_param == "force")
//     {
//       user_control.force_override = true;
//     }
//     else if (current_param == "fi")
//     {
//       user_control.initialiser_path = current_value;
//     }
//     break;
//   }
//   case 's':
//   {

//     if (current_param == "serde")
//     {
//       user_control.serde = true;
//       user_control.serde_file = std::string(current_value);
//     }
//     break;
//   }

//   default:
//   {
//     throw_bad_arg(current_param);
//     break;
//   }
//   }
// }

static CliResults<Core::UserControlParameters>
sanitise_check_cli(Core::UserControlParameters&& params)
{

  if (params.delta_time < 0)
  {
    return CliResults<Core::UserControlParameters>("Wrongtime step (d_t<0)");
  }

  if (params.number_particle == 0)
  {
    return CliResults<Core::UserControlParameters>(
        "Missing number of particles");
  }

  if (params.final_time <= 0)
  {
    return CliResults<Core::UserControlParameters>(
        "Final time must be positive");
  }

  return CliResults<Core::UserControlParameters>(std::move(params));
}

std::string get_help_message() noexcept
{
  std::stringstream os;
  os << "Usage: ";
  print_red(os, "BIOCMA-MCST");
  // os << "  -np <number_of_particles> [-ff <flow_file_folder_path>] [OPTIONS]
  // " << '\n'; os << "\nMandatory arguments:" << '\n'; os << "  -np <number>,
  // --number-particles <number>\tNumber of particles" << '\n'; os << "  -d
  // <number>, --duration <number>\tSimulation duration" << '\n'; os << "  -ff
  // <flow_file_folder_path>\t\tPath to flow file (Default: "
  //       "./rawdata)"
  //    << '\n';

  // os << "\nOptional arguments:" << '\n';
  // os << "  -h, --help\t\tDisplay this help message" << '\n';
  // os << "  -v, --verbose\t\tVerbose mode" << '\n';
  // os << "  -nt <number>, --number-threads <number>\tNumber of threads per "
  //       "process"
  //    << '\n';

  // os << "Available model:\r\n";
  // os << "WIP\r\n";
  // for (const auto& i : get_available_models())
  // {
  //   os << i << "\r\n";
  // }

  os << "\nExample:" << '\n';
  os << "  BIOCMA-MCST -np 100 -ff /path/to/flow_file_folder/ [-v]" << '\n';

  return os.str();
}

// static void print_green(std::ostream &os, std::string_view message)
// {
//   os << "\033[1;32m" << message << "\033[0m";
// }

static void print_red(std::ostream& os, std::string_view message)
{
  os << "\033[1;31m" << message << "\033[0m"; // ANSI escape code for red color
}

static void throw_bad_arg(std::string_view arg)
{
  std::string msg = "Unknown argument: ";
  msg += std::string(arg);
  throw std::invalid_argument(msg);
}
