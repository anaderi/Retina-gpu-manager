#include "LocalSocketClient.h"

#include "EofException.h"
#include "SystemException.h"

#include <cstring>
#include <stdexcept>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

#include <sstream>

//-------------
// construction
//-------------

LocalSocketClient::LocalSocketClient(const char * path) :
    m_path   (path),
    m_socket (-1) {
  m_socket = socket(AF_UNIX, SOCK_STREAM, 0);
  if (m_socket == -1)
    throw SystemException("Could not create socket.");

  sockaddr_un address = { 0, { 0 } };
  address.sun_family = AF_LOCAL;
  const size_t maxPathLength = sizeof(address.sun_path) - 1;
  if (m_path.size() > maxPathLength)
    throw std::runtime_error("The socket path is too long.");
  std::strncpy(address.sun_path, m_path.c_str(), maxPathLength);

  if (-1 == connect(m_socket, (sockaddr*)&address, SUN_LEN(&address)))
    throw SystemException("Could not connect to socket.");
}

LocalSocketClient::~LocalSocketClient() {
  if (m_socket != -1)
    close(m_socket);
}

//--------------------------
// ITransport implementation
//--------------------------

void LocalSocketClient::readBytes(void * data, size_t size) {
  size_t total = 0u;
  while (total < size) {
    ssize_t received = read(m_socket, static_cast<uint8_t*>(data) + total, size - total);
    if (received == 0)
      throw EofException();
    if (received < 0)
      throw SystemException("Read error.");
    total += static_cast<size_t>(received);
  }
}

void LocalSocketClient::writeBytes(const void * data, size_t size) {
  size_t total = 0u;
  while (total < size) {
    ssize_t sent = write(m_socket, (uint8_t*)data + total, size - total);
    if (sent < 0)
      throw SystemException("Write error.");
    total += static_cast<size_t>(sent);
  }
}
