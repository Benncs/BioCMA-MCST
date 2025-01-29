#ifndef __TOOL_EMBED_HPP__
#define __TOOL_EMBED_HPP__
#include <span>
#include <string_view>

namespace EmbedPython
{
  void call_add_function(std::span<double> data, size_t n_r, size_t n_c,std::span<double> datac, size_t cn_r, size_t cn_c);
  void save_pickle(std::string_view name, std::span<double> data, size_t n_r, size_t n_c);
  void call_show_function();
  void execute_script(const std::string& script);

  void execute_script_from_file(const std::string& filename);

  

} // namespace EmbedPython

#endif