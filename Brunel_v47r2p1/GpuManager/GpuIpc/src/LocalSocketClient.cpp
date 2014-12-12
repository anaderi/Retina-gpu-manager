#include "LocalSocketClient.h"
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
    size_t received = read(m_socket, static_cast<uint8_t*>(data) + total, size - total);
    if (received == static_cast<size_t>(-1)) {
      std::stringstream msg;
      msg << "Read error; " << total << " bytes received.";
      throw std::runtime_error(msg.str().c_str());
    }
    total += received;
  }
  if (total > size) {
    std::stringstream msg;
    msg << "Asked for " << size << " bytes; received " << total << ".";
    throw std::runtime_error(msg.str().c_str());
  }
}

void LocalSocketClient::writeBytes(const void * data, size_t size) {
  size_t sent = write(m_socket, data, size);
  if (sent == static_cast<size_t>(-1))
    throw SystemException("Write error.");
  if (sent != size) {
    std::stringstream msg;
    msg << "Asked for " << size << " bytes; sent " << sent << ".";
    throw std::runtime_error(msg.str().c_str());
  }
}
