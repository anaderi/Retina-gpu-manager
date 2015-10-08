#pragma once

#include <string>

#include "IConnector.h"

/// Connector for the SocketServer.
class LocalSocketServerConnector : public IConnector {
  public: // public interface

    /// Create a connector for a socket with the given path.
    LocalSocketServerConnector(const char * path);

    virtual ~LocalSocketServerConnector();

  public: // IConnector implementation

    /// Block until a client connects or throw an exception.
    virtual std::shared_ptr<ITransport> accept();

    /// Called when the server is halted.
    virtual void close();

  private: // data

    std::string m_path;
    int         m_socket;
};
