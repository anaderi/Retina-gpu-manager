/*
 * Tcp socket server adapted from the LocalSocketServerConnector
 * by Alexey Badalov and Marco Corvo
 */

#include "SocketServer.h"
#include "TcpSocketServerConnector.h"
#include "SystemException.h"

#include <arpa/inet.h>
#include <cstring>
#include <netinet/in.h>
#include <stdexcept>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

//-----------------
// public interface
//-----------------

TcpSocketServerConnector::TcpSocketServerConnector(int port, const std::string& hname) :
    m_host   (hname),
    m_port   (port),
    m_socket (-1) {
  m_socket = ::socket(AF_INET, SOCK_STREAM, 0);
  if (m_socket == -1)
    throw SystemException("Could not create socket.");

  sockaddr_in serv_addr = {};

  serv_addr.sin_family = AF_INET;
  if (m_host.compare("localhost") == 0)
      serv_addr.sin_addr.s_addr = INADDR_ANY;
  else
      serv_addr.sin_addr.s_addr = inet_addr(m_host.c_str());
  serv_addr.sin_port = htons(m_port);
  if (::bind(m_socket, (sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
      throw SystemException("Could not bind to socket");
  }
  if (::listen(m_socket,5) < 0) {
      throw SystemException("Could not listen to socket");
  }

}

TcpSocketServerConnector::~TcpSocketServerConnector() {
  if (m_socket != -1)
    ::close(m_socket);
}

//--------------------------
// IConnector implementation
//--------------------------

std::shared_ptr<ITransport> TcpSocketServerConnector::accept() {
  int accepted = ::accept(m_socket, NULL, NULL);
  if (-1 == accepted)
    throw SystemException("Could not accept incoming connection.");
  return std::make_shared<SocketServer>(accepted);
}

void TcpSocketServerConnector::close() {
  shutdown(m_socket, SHUT_RDWR);
}
