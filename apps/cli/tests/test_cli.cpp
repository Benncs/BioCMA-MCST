#include <cassert>
#include <cli_parser.hpp>

int main()
{
  {
    const char* argv[] = {"program_name",
                          "-er",
                          "result.txt",
                          "-mn",
                          "model_1",
                          "-nt",
                          "8",
                          "-np",
                          "1000000",
                          "-d",
                          "3.14",
                          "-r"};
    int argc = 9;
    auto result = parse_cli(argc, const_cast<char**>(argv));

    //     assert(result.has_value()); // Should return a valid UserControlParameters object
    //     const auto& control = result.value();

    //     // Validate values are correctly set
    //     assert(control.results_file_name == "result.txt");
    //     assert(control.model_name == "model_1");
    //     assert(control.n_thread == 8);
    //     assert(control.number_particle == 1000000);
    //     assert(control.final_time == 3.14);
    //     assert(control.recursive == true);
    //   }

    //   {
    //     const char* argv[] = {"program_name"};
    //     int argc = 1;
    //     auto result = parse_cli(argc, const_cast<char**>(argv));

    //     assert(result == std::nullopt); // Should return std::nullopt because of insufficient
    //     arguments
  }
}