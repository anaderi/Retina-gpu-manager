#pragma once

#include "ITransport.h"

#include <string>

/// Transport over a socket.
class SocketServer : public ITransport {
public: // public interface

  /// Construct a server instance for a given socket handle.
  SocketServer(int socket);

  virtual ~SocketServer();

public: // ITransport implementation

  /// Read a block of bytes.
  virtual void readBytes(void * data, size_t size);

  /// Write a block of bytes.
  virtual void writeBytes(const void * data, size_t size);

private: // data

  int m_socket;
};
