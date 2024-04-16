#ifndef __CLI_PARSER_HPP__
#define __CLI_PARSER_HPP__

#include <common/common.hpp>
#include <iostream>

SimulationParameters parse_cli(int argc, char **argv);
void showHelp(std::ostream &os);

#endif //__CLI_PARSER_HPP__
