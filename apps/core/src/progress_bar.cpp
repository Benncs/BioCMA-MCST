#include <iomanip>
#include <iostream>
#include <progress_bar.hpp>

namespace IO
{
  constexpr size_t PROGRESS_BAR_WIDTH = 100; // Number of
  constexpr char PROGRESS_BAR_SYMBOL = '*';

  ProgressBar::ProgressBar()
      : default_precision(
            std::cout.precision()) // NOLINT Conversion long to int)
  {
    buffer = std::string(PROGRESS_BAR_WIDTH, ' ');
  }

  void ProgressBar::show(std::ostream& out_stream,
                         size_t total,
                         size_t current_position)
  {

    std::ios::sync_with_stdio(false);
    size_t progress = (current_position * PROGRESS_BAR_WIDTH) / total;

    std::fill_n(buffer.begin(), progress, PROGRESS_BAR_SYMBOL);
    out_stream << "Progress: [" << buffer << "] " << std::fixed
               << std::setprecision(2)
               << (static_cast<float>(current_position) * 100.0 /
                   static_cast<float>(total))
               << "%\r" << std::flush << std::setprecision(default_precision);
  }

  void ProgressBar::show_percentage(std::ostream& out_stream,
                                    size_t total,
                                    size_t current_position) const
  {
    std::ios::sync_with_stdio(false);
    out_stream << "Progress: [" << std::fixed << std::setprecision(2)
               << (static_cast<float>(current_position) * 100.0 /
                   static_cast<float>(total))
               << "%" << "]\r" << std::flush
               << std::setprecision(default_precision);
  }

} // namespace IO
