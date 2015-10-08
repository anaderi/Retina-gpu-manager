#include "TcpSocketClient.h"

#include "EofException.h"
#include "SystemException.h"

#include <arpa/inet.h>
#include <cstring>
#include <netinet/in.h>
#include <stdexcept>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

#include <sstream>

//-------------
// construction
//-------------

TcpSocketClient::TcpSocketClient(int port, const std::string& hname) :
    m_host   (hname),
    m_socket (-1),
    m_port   (port) {
  m_socket = socket(AF_INET, SOCK_STREAM, 0);
  if (m_socket == -1)
    throw SystemException("Could not create socket.");

  sockaddr_in address = {};
  address.sin_family      = AF_INET;
  address.sin_port        = htons(m_port);
  address.sin_addr.s_addr = inet_addr(m_host.c_str());

  if (-1 == connect(m_socket, (sockaddr*)&address, sizeof(address)))
    throw SystemException("Could not connect to socket.");
}

TcpSocketClient::~TcpSocketClient() {
  if (m_socket != -1)
    close(m_socket);
}

//--------------------------
// ITransport implementation
//--------------------------

void TcpSocketClient::readBytes(void * data, size_t size) {
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

void TcpSocketClient::writeBytes(const void * data, size_t size) {
  size_t total = 0u;
  while (total < size) {
    ssize_t sent = write(m_socket, (uint8_t*)data + total, size - total);
    if (sent < 0)
      throw SystemException("Write error.");
    total += static_cast<size_t>(sent);
  }
}
