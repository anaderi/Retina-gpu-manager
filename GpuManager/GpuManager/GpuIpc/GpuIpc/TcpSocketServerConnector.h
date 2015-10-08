#pragma once

#include <string>

#include "IConnector.h"

/// Connector for the SocketServer.
class TcpSocketServerConnector : public IConnector {
  public: // public interface

    /// Create a connector for a socket with the given path.
    TcpSocketServerConnector(int port, const std::string& hname);

    virtual ~TcpSocketServerConnector();

  public: // IConnector implementation

    /// Block until a client connects or throw an exception.
    virtual std::shared_ptr<ITransport> accept();

    /// Called when the server is halted.
    virtual void close();

  private: // data

    const std::string m_host;
    int   m_port;
    int   m_socket;
};
