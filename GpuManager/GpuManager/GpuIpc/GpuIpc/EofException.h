#pragma once

#include <stdexcept>

/// Thrown when the other side hangs up during communication.
class EofException : public std::runtime_error {
public:
  EofException() : runtime_error("End of file.") {}
  virtual ~EofException() noexcept {}
};
