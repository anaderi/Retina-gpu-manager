#pragma once

#include <stdexcept>
#include <string>

/// Serious exception, recording errno.
class SystemException : public std::runtime_error {
public:

  /// Construct an instance of SystemException with the given message.
  explicit SystemException(const std::string & message);

  /// Construct an instance of SystemException with the given message and error ID.
  SystemException(const std::string & message, int error);

  virtual ~SystemException() noexcept;

  /// Return the errno corresponding to the error.
  int error() const noexcept;

private:

  int m_error;
};
