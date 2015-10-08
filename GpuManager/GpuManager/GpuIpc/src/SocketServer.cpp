#include "SocketServer.h"

#include "EofException.h"
#include "SystemException.h"

#include <cstring>
#include <stdint.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

//-----------------
// public interface
//-----------------

SocketServer::SocketServer(int socket) :
    m_socket (socket) {
}

SocketServer::~SocketServer() {
  if (m_socket != -1)
    close(m_socket);
}

//--------------------------
// ITransport implementation
//--------------------------

void SocketServer::readBytes(void * data, size_t size) {
  size_t total = 0u;
  while (total < size) {
    ssize_t received = read(m_socket, (uint8_t*)data + total, size - total);
    if (received == 0)
      throw EofException();
    if (received < 0)
      throw SystemException("Read error.");
    total += static_cast<size_t>(received);
  }
}

void SocketServer::writeBytes(const void * data, size_t size) {
  size_t total = 0u;
  while (total < size) {
    ssize_t sent = write(m_socket, (uint8_t*)data + total, size - total);
    if (sent < 0)
      throw SystemException("Write error.");
    total += static_cast<size_t>(sent);
  }
}
