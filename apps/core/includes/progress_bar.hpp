#ifndef __IO_PROGRESS_BAR_HPP__
#define __IO_PROGRESS_BAR_HPP__
#include <ostream>

namespace IO
{
  class ProgressBar
  {
  public:
    ProgressBar();
    void show(std::ostream& out_stream,
              std::size_t total,
              std::size_t current_position);

    void show_percentage(std::ostream& out_stream,
                         std::size_t total,
                         std::size_t current_position) const;

  private:
    std::string buffer;

    const int default_precision; // NOLINT
  };
} // namespace IO

#endif
