#include <cli_parser.hpp>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string_view>
// static void parseOptional(SimulationParameters& params, std::string_view
// arg);
static void check_cli(SimulationParameters &params);
// static void print_green(std::ostream &os, std::string_view message);
static void print_red(std::ostream &os, std::string_view message);
static void throw_bad_arg(std::string_view arg);

static void parseArg(SimulationParameters &params,
                     std::string current_param,
                     std::string_view current_value,
                     bool &recur);

static void recur_path(std::string_view rootPath, SimulationParameters &params)
{
  size_t count = 1;
  std::string dirName = "i_" + std::to_string(count) + "/";
  std::filesystem::path dirPath = std::string(rootPath) + dirName;
  while (std::filesystem::exists(dirPath) &&
         std::filesystem::is_directory(dirPath))
  {
    ++count;
    params.flow_files.push_back(dirPath.string());
    dirName = "i_" + std::to_string(count) + "/";
    dirPath = std::string(rootPath) + dirName;
  }
}

std::optional<SimulationParameters> parse_cli(int argc, char **argv) noexcept
{
  SimulationParameters params = SimulationParameters::m_default();

  params.flow_files.clear();

  bool recur = false;
  try
  {
    int iarg = 1;
    while (iarg < argc - 1)
    {
      auto current_param = std::string(argv[iarg]);
      auto current_value = std::string_view(argv[iarg + 1]);
      if (current_param.data() != nullptr && current_param[0] != '\0')
      {

        parseArg(params, current_param, current_value, recur);
      }
      else
      {
        throw_bad_arg("");
      }
      iarg += 2;
    }
    auto root = params.flow_files[0];
         params.root = std::string(root);
    if (params.flow_files.size() == 1 && recur)
    {
      auto root = params.flow_files[0];
      params.flow_files.clear();
      recur_path(root, params);
    }
    check_cli(params);
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

  return params;
}

static void parseArg(SimulationParameters &params,
                     std::string current_param,
                     std::string_view current_value,
                     bool &recur)
{
  std::string path;
  current_param = std::string(current_param.begin() + 1, current_param.end());
  switch (current_param[0])
  {
  case 'n':
  {
    if (current_param == "nt")
    {
      params.n_threads = std::stoi(std::string(current_value));
    }
    else if (current_param == "np")
    {
      params.n_particles = std::stol(std::string(current_value));
    }
    break;
  }
  case 'd':
  {
    if (current_param == "dt")
    {
      params.d_t = std::stod(std::string(current_value));
    }
    else if (current_param == "d")
    {
      params.final_time = std::stod(std::string(current_value));
    }
    break;
  }
  case 'r':
  {
    recur = true;
    break;
  }
  case 'f':
  {
    params.flow_files.emplace_back(current_value);

    break;
  }

  default:
    throw_bad_arg(current_param);
    break;
  }
}

// static void parseOptional(SimulationParameters& params, std::string_view arg)
// {
//     arg = std::string_view(arg.begin() + 1, arg.end()); // Resize view to
//     skip '-' if (arg.size() == 0) {
//         throw_bad_arg(""); // Throw an exception with an appropriate error
//         message
//     }

//     switch (arg[0]) {
//     case '-': {
//         // Fullname options
//         if (arg == "-help") {
//             showHelp(std::cout);
//             exit(0);
//         } else if (arg == "-init") {
//         } else if (arg == "-duration") {
//             params.final_time = std::stod(/*Split space*/);
//         } else if (arg == "-number-threads") {
//             params.n_threads = std::stoi(/*Split space*/);
//         } else if (arg == "-number-particles") {
//             params.n_particles = std::stol(/*Split space*/);
//         } else {
//             throw_bad_arg(arg);
//         }
//         break;
//     }
//     case 'n': {
//         if (arg.size() == 2) {
//             if (arg[1] == 'p') {
//                 // Handle number of particles option
//                 params.n_particles = std::stol(std::string(arg.begin() + 2,
//                 arg.end()));
//             } else {
//                 throw_bad_arg(arg); // Throw an exception with the unknown
//                 option
//             }
//         } else {
//             throw_bad_arg(arg); // Throw an exception with the unknown option
//         }
//         break;
//     }
//     case 't': {
//         // Handle number of threads option
//         params.n_threads = std::stoi(std::string(arg.begin() + 1,
//         arg.end())); break;
//     }
//     case 'i': {
//         // Handle 'i' option if needed
//         break;
//     }
//     case 'h': {
//         showHelp(std::cout);
//         exit(0);
//     }
//     default: {
//         throw_bad_arg(arg);
//         break;
//     }
//     }
// }

static void check_cli(SimulationParameters &params)
{
  if (params.flow_files.empty())
  {
    throw std::invalid_argument("Missing files path");
  }

  if (params.d_t < 0)
  {
    throw std::invalid_argument("Wrongtime step (d_t<0)");
  }

  if (params.n_particles == 0)
  {
    throw std::invalid_argument("Missing number of particles");
  }

  if (params.final_time <= 0)
  {
    throw std::invalid_argument("Final time must be positive");
  }
}

void showHelp(std::ostream &os)
{
  os << "Usage: ";
  print_red(os, "BIOCMA-MCST");
  os << "  -np <number_of_particles> [-ff <flow_file_folder_path>] [OPTIONS] "
     << '\n';
  os << "\nMandatory arguments:" << '\n';
  os << "  -np <number>, --number-particles <number>\tNumber of particles"
     << '\n';
  os << "  -d <number>, --duration <number>\tSimulation duration" << '\n';
  os << "  -ff <flow_file_folder_path>\t\tPath to flow file (Default: "
        "./rawdata)"
     << '\n';

  os << "\nOptional arguments:" << '\n';
  os << "  -h, --help\t\tDisplay this help message" << '\n';
  os << "  -v, --verbose\t\tVerbose mode" << '\n';
  os << "  -nt <number>, --number-threads <number>\tNumber of threads per "
        "process"
     << '\n';

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
  msg += arg;
  throw std::invalid_argument(msg);
}
