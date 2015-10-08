#include "AdminClient.h"
#include "AdminID.h"
#include "GpuIpc/IProtocol.h"

using namespace std;

AdminClient::AdminClient(IProtocol & protocol) : m_protocol(protocol) {
}

void AdminClient::exit() {
  m_protocol.writeUInt32(AdminID_Exit);
}

void AdminClient::loadHandler(const string & handlerName) {
  m_protocol.writeUInt32(AdminID_LoadHandler);
  m_protocol.writeString(handlerName);
}
