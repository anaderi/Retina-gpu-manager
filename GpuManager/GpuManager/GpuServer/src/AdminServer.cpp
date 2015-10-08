#include "AdminServer.h"
#include "AdminID.h"
#include "GpuIpc/IProtocol.h"
#include "IApp.h"

#include <stdexcept>

AdminServer::AdminServer(IApp & app) : m_app (app) {
}

void AdminServer::process(IProtocol & protocol) {
  uint32_t id = protocol.readUInt32();
  switch (id) {
  case AdminID_Exit:        processExit(protocol);        break;
  case AdminID_LoadHandler: processLoadHandler(protocol); break;
  default:
    throw std::runtime_error("AdminServer::process: unknown ID");
  }
}

void AdminServer::processExit(IProtocol &) {
    m_app.exit();
}

void AdminServer::processLoadHandler(IProtocol & protocol) {
  m_app.loadHandler(protocol.readString());
}
