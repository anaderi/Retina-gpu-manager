#pragma once

#include "ITransport.h"

#include <string>

/// Transport over a named socket.
class TcpSocketClient : public ITransport {
  public:
    /// Connect to the socket with the given path.
    TcpSocketClient(int port, const std::string& hname);

    virtual ~TcpSocketClient();

  public: // ITransport implementation

    /// Read a block of bytes.
    virtual void readBytes(void * data, size_t size);

    /// Write a block of bytes.
    virtual void writeBytes(const void * data, size_t size);

  private:
    std::string m_host;
    int         m_socket;
    int         m_port;
};
