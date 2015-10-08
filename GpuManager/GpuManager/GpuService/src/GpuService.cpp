#include "GpuService.h"

#include "GpuIpc/LocalSocketClient.h"
#include "GpuIpc/Protocol.h"
#include "GpuIpc/TcpSocketClient.h"

#include <GaudiKernel/SvcFactory.h>

#include <stdexcept>
#include <string>

using namespace std;

DECLARE_SERVICE_FACTORY(GpuService)

//-------------
// construction
//-------------

GpuService::GpuService(const std::string & name, ISvcLocator * sl) :
    base_class        (name, sl),
    m_transport       (nullptr),
    m_protocol        (nullptr),
    m_connectionType  ("local"),
    m_localSocketPath ("/tmp/GpuManager"),
    m_serverHost      ("localhost"),
    m_serverPort      (65000) {
  declareProperty("ConnectionType", m_connectionType);
  declareProperty("SocketPath",     m_localSocketPath);
  declareProperty("TcpHost",        m_serverHost);
  declareProperty("TcpPort",        m_serverPort);
}

GpuService::~GpuService() {
}

//-----------------------------
// IGpuService implementation
//-----------------------------

void GpuService::submitData(
    std::string  handlerName,
    const void * data,
    const size_t size,
    Alloc        allocResults,
    AllocParam   allocResultsParam) {
  // send the name of the addressee
  m_protocol->writeString(handlerName);

  // send the data package
  m_protocol->writeUInt32(size);
  m_protocol->writeData(data, size);

  // receive the result
  size_t resultSize = m_protocol->readUInt32();

  // handle errors
  const size_t FAIL_FLAG = 0xFFFFFFFF;
  if (resultSize == FAIL_FLAG) {
    string message = m_protocol->readString();
    throw runtime_error(message);
  }

  // zero data size is normal
  if (resultSize > 0) {
    void * resultData = allocResults(resultSize, allocResultsParam);
    if (resultData == NULL) {
      // allocate a temporary buffer, if the client won't cooperate
      vector<uint8_t> temp(resultSize);
      m_protocol->readData(&temp[0], resultSize);
    } else {
      m_protocol->readData(resultData, resultSize);
    }
  }

  // receive and discard performance data
  m_protocol->readDouble();
}

//-----------------------
// Service implementation
//-----------------------

StatusCode GpuService::initialize() {
  StatusCode sc(Service::initialize());
  if (StatusCode::SUCCESS != sc)
    return sc;

  acquireIO();

  return StatusCode::SUCCESS;
}

StatusCode GpuService::finalize() {
  releaseIO();
  return Service::finalize();
}

//------------------
// private functions
//------------------

void GpuService::acquireIO() {
  if (m_transport)
    return;
  if (m_connectionType.value() == "local") {
    m_transport = new LocalSocketClient(m_localSocketPath.value().c_str());
    m_protocol  = new Protocol(*m_transport);
  } else if (m_connectionType.value() == "tcp") {
    m_transport = new TcpSocketClient(m_serverPort, m_serverHost);
    m_protocol  = new Protocol(*m_transport);
  } else {
    throw runtime_error("Invalid connection type. Should be one of: local, tcp.");
  }
}

void GpuService::releaseIO() {
  if (m_protocol != nullptr) {
    delete m_protocol;
    m_protocol = nullptr;
  }
  if (m_transport != nullptr) {
    delete m_transport;
    m_transport = nullptr;
  }
}
