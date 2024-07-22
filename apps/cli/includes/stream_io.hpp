#ifndef __STREAM_IO_HPP__   
#define __STREAM_IO_HPP__


#include <streambuf>


int redirect_stdout(std::streambuf *&original_buffer,
                           std::stringstream &variable_stream);
void restore_stdout(int original_stdout,
                           std::streambuf *&original_buffer);


#define REDIRECT_BLOCK(__block__, __f_verbose__, __f__redirect)                \
  if constexpr ((__f__redirect))                                               \
  {                                                                            \
    std::stringstream output_variable;                                         \
    output_variable.str("");                                                   \
    std::streambuf *original_buffer = nullptr;                                 \
    auto original_stdout_fd =                                                  \
        redirect_stdout(original_buffer, output_variable);                     \
    {__block__} restore_stdout(original_stdout_fd, original_buffer);           \
    if constexpr ((__f_verbose__))                                             \
    {                                                                          \
      std::cout << output_variable.str() << std::endl;                         \
    }                                                                          \
  }                                                                            \
  else                                                                         \
  {                                                                            \
    __block__                                                                  \
  }

#endif 