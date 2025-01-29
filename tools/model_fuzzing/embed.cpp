#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <span>
#include <string_view>
#include <fstream>
#include <iostream> 

namespace EmbedPython
{

  namespace py = pybind11;

  void save_pickle(std::string_view name, std::span<double> data, size_t n_r, size_t n_c)
  {
    if (data.size() != n_r * n_c)
    {
      throw std::runtime_error("Data size does not match specified dimensions.");
    }

    py::array_t<double> arr({n_r, n_c}, data.data());
    py::module pickle = py::module::import("pickle");
    py::bytes pickled_data = pickle.attr("dumps")(arr);
    std::ofstream out_file(name.data(), std::ios::binary);
    out_file.write(pickled_data.cast<std::string>().c_str(),
                   pickled_data.cast<std::string>().size());
    out_file.close();
  }

  void call_add_function(std::span<double> data, size_t n_r, size_t n_c,std::span<double> datac, size_t cn_r, size_t cn_c)
  {
    py::array_t<double> arr({n_r, n_c}, data.data());
    
    py::array_t<double> c({cn_r, cn_c}, datac.data());


    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("tools/model_fuzzing");
    py::module show_module = py::module::import("fuzzing_plot");
    show_module.attr("add_fig")(arr,c);
  }

   void call_show_function()
  {
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("tools/model_fuzzing");
    py::module show_module = py::module::import("fuzzing_plot");
    show_module.attr("show")();
  }

  void execute_script(const std::string& script)
  {
    try
    {
      // Execute the script
      py::exec(script);
      std::cout << "Python script executed successfully." << std::endl;
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error executing script: " << e.what() << std::endl;
    }
  }

  void execute_script_from_file(const std::string& filename)
  {
    // Read the Python script from a file
    std::ifstream script_file(filename);
    if (!script_file.is_open())
    {
      std::cerr << "Failed to open script file: " << filename << std::endl;
      return;
    }

    std::string script((std::istreambuf_iterator<char>(script_file)),
                       std::istreambuf_iterator<char>());

    // Execute the script
    execute_script(script);
  }

} // namespace EmbedPython


