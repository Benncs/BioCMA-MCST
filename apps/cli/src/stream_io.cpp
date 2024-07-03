#include <iostream>
#include <stream_io.hpp>
#include <unistd.h>

#include <sstream>

static bool is_stdout_redirect = false;

void restore_stdout(int original_stdout_fd,
                           std::streambuf *&original_buffer)
{
  if (is_stdout_redirect)
  {
    std::cout.rdbuf(original_buffer);
    fflush(stdout); // Flush the buffer to ensure all previous output is written
    dup2(original_stdout_fd,
         fileno(stdout)); // Restore the original file descriptor of stdout

    is_stdout_redirect = false;
  }
}

int redirect_stdout(std::streambuf *&original_buffer,
                           std::stringstream &variable_stream)
{
  // Check if redirection is already active
  if (!is_stdout_redirect)
  {
    // Save the original buffer of std::cout
    original_buffer = std::cout.rdbuf();
    // Redirect std::cout to the stringstream
    std::cout.rdbuf(variable_stream.rdbuf());

    fflush(stdout); // Flush the buffer to ensure all previous output is written
    int original_stdout_fd =
        dup(fileno(stdout)); // Save the original file descriptor of stdout
    auto *fd = freopen("/dev/null", "w", stdout); // Redirect stdout to /dev/null

    // Set the flag to indicate that redirection is active
    is_stdout_redirect = true;
    return original_stdout_fd;
  }
  return fileno(stdout);
}


                           