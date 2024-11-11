#include "rt_init.hpp"
#include <cli_parser.hpp>
#include <core/simulation_parameters.hpp>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string_view>

static void sanitise_check_cli(const Core::UserControlParameters &params);
// static void print_green(std::ostream &os, std::string_view message);
static void print_red(std::ostream &os, std::string_view message);
static void throw_bad_arg(std::string_view arg);

static void
parseArg(Core::UserControlParameters &user_controll, std::string current_param, std::string_view current_value);


static std::optional<Core::UserControlParameters> parse_user_param(int argc, char **argv)
{
  Core::UserControlParameters control = Core::UserControlParameters::m_default();

  try
  {
    int iarg = 1;
    while (iarg < argc - 1)
    {
      auto current_param = std::string(argv[iarg]);
      auto current_value = std::string_view(argv[iarg + 1]);
      if (current_param.data() != nullptr && current_param[0] != '\0')
      {
        if (current_param == "h")
        {
          showHelp(std::cout);
          exit(0);
        }

        parseArg(control, current_param, current_value);
      }
      else
      {
        throw_bad_arg("ALLO");
      }
      iarg += 2;
    }
    return control;
  }
  catch (std::invalid_argument &e)
  {
    std::cerr << e.what() << '\n';
    return std::nullopt;
  }
  catch (std::exception const &e)
  {
    std::cerr << e.what() << '\n';
    return std::nullopt;
  }
}

std::optional<Core::UserControlParameters> parse_cli(int argc, char **argv) noexcept
{
  auto opt_control = parse_user_param(argc, argv);
  if (!opt_control.has_value())
  {
    return std::nullopt;
  }
  auto control = *opt_control;
  sanitise_check_cli(control);

  return control;
}

static void
parseArg(Core::UserControlParameters &user_control, std::string current_param, std::string_view current_value)
{
  std::string path;
  // TODO need check that begin()+1 is not OOB
  current_param = std::string(current_param.begin() + 1, current_param.end());
  switch (current_param[0])
  {
  case 'e':
  {
    if (current_param == "er")
    {
      user_control.results_file_name = std::string(current_value);
    }
    break;
  }
  case 'm':
  {
    if (current_param == "mn")
    {
      user_control.model_name = std::string(current_value);
    }
    break;
  }
  case 'n':
  {
    if (current_param == "nt")
    {
      user_control.n_thread = std::stoi(std::string(current_value));
    }
    else if (current_param == "np")
    {
      user_control.number_particle = std::stol(std::string(current_value));
    }
    else if (current_param == "nex")
    {
      user_control.number_exported_result = std::stol(std::string(current_value));
    }
    break;
  }
  case 'd':
  {
    if (current_param == "dt")
    {
      user_control.delta_time = std::stod(std::string(current_value));
    }
    else if (current_param == "d")
    {
      user_control.final_time = std::stod(std::string(current_value));
    }
    break;
  }

  case 'r':
  {
    user_control.recursive = true;
    break;
  }
  case 'f':
  {
    if (current_param == "f")
    {
      user_control.cma_case_path = current_value;
    }
    else if (current_param == "force")
    {
      user_control.force_override = true;
    }
    else if (current_param == "fi")
    {
      user_control.initialiser_path = current_value;
    }
    break;
  }
  case 's':
  {

    if (current_param == "serde")
    {
      user_control.serde = true;
      user_control.serde_file = std::string(current_value);
    }
    break;
  }

  default:
  {
    throw_bad_arg(current_param);
    break;
  }
  }
}

static void sanitise_check_cli(const Core::UserControlParameters &params)
{

  if (params.delta_time < 0)
  {
    throw std::invalid_argument("Wrongtime step (d_t<0)");
  }

  if (params.number_particle == 0)
  {
    throw std::invalid_argument("Missing number of particles");
  }

  if (params.final_time <= 0)
  {
    throw std::invalid_argument("Final time must be positive");
  }
}

void showHelp(std::ostream &os) noexcept
{
  os << "Usage: ";
  print_red(os, "BIOCMA-MCST");
  // os << "  -np <number_of_particles> [-ff <flow_file_folder_path>] [OPTIONS] " << '\n';
  // os << "\nMandatory arguments:" << '\n';
  // os << "  -np <number>, --number-particles <number>\tNumber of particles" << '\n';
  // os << "  -d <number>, --duration <number>\tSimulation duration" << '\n';
  // os << "  -ff <flow_file_folder_path>\t\tPath to flow file (Default: "
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
}

// static void print_green(std::ostream &os, std::string_view message)
// {
//   os << "\033[1;32m" << message << "\033[0m";
// }

static void print_red(std::ostream &os, std::string_view message)
{
  os << "\033[1;31m" << message << "\033[0m"; // ANSI escape code for red color
}

static void throw_bad_arg(std::string_view arg)
{
  std::string msg = "Unknown argument: ";
  msg += std::string(arg);
  throw std::invalid_argument(msg);
}
