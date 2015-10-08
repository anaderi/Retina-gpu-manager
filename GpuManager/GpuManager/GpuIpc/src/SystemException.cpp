#include "SystemException.h"

#include <cerrno>
#include <cstring>
#include <sstream>

using namespace std;

namespace {
  string makeMessage(int error, const string & message) {
    stringstream msg;
    msg << message << " (" << error << ") " << strerror(error);
    return msg.str();
  }
}

SystemException::SystemException(const string & message) :
    SystemException(message, errno) {
}

SystemException::SystemException(const string & message, int error) :
    runtime_error (makeMessage(error, message)),
    m_error       (error) {
}

SystemException::~SystemException() noexcept {
}

int SystemException::error() const noexcept {
  return m_error;
}

