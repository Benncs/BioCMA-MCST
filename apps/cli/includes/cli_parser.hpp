#ifndef __CLI_PARSER_HPP__
#define __CLI_PARSER_HPP__

#include <common/common.hpp>
#include <iostream>
#include <optional>

std::optional<SimulationParameters> parse_cli(int argc, char **argv) noexcept;
void showHelp(std::ostream &os);

#endif //__CLI_PARSER_HPP__
