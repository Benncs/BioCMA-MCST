#include <cli_parser.hpp>
#include <iostream>

SimulationParameters parseCLI(int argc, char **argv)
{

  std::cout<<argc<<std::endl;
  auto file = "/home/benjamin/Documenti/code/cpp/BIREM_generate/out/";
  return {1'000, 3, 10., {file}};
}