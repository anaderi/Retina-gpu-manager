#pragma once

#include "ITransport.h"

#include <string>

/// Transport over a named socket.
class LocalSocketClient : public ITransport {
  public:
    /// Connect to the socket with the given path.
    LocalSocketClient(const char * path);

    virtual ~LocalSocketClient();

  public: // ITransport implementation

    /// Read a block of bytes.
    virtual void readBytes(void * data, size_t size);

    /// Write a block of bytes.
    virtual void writeBytes(const void * data, size_t size);

  private:
    std::string m_path;
    int         m_socket;
};
