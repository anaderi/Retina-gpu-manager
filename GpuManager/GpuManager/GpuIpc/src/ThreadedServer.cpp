#include "ThreadedServer.h"

#include "EofException.h"
#include "IConnector.h"
#include "IProcessor.h"
#include "IProtocol.h"
#include "ITransport.h"
#include "SystemException.h"
#include "ThreadGroup.h"

#include <syslog.h>

#include <iostream>
#include <utility>

using namespace std;

namespace {
  void printError(const char * msg) {
      cerr << msg << '\n';
      syslog(LOG_ERR, "%s", msg);
  }
}

//-----------------
// public interface
//-----------------

ThreadedServer::ThreadedServer(
    IConnector      & connector,
    ProtocolFactory   protocolFactory,
    IProcessor      & processor) :
    // initializers
    m_threadID        (0),
    m_connector       (connector),
    m_protocolFactory (protocolFactory),
    m_processor       (processor) {
}

boost::thread * ThreadedServer::serve() {
  return new boost::thread(&ThreadedServer::serveConnector, this);
}

void ThreadedServer::stop() {
  m_connector.close();
}

//------------------
// private functions
//------------------

void ThreadedServer::serveConnector() {
  ThreadGroup connections;
  try {
    for (;;) { // until an exception is thrown
      std::shared_ptr<ITransport> transport(m_connector.accept());
      std::shared_ptr<IProtocol>  protocol(m_protocolFactory(*transport));
      connections.createThread(ConnectionHandler(*this, transport, protocol, m_processor));
    }
  } catch (const SystemException & e) {
    switch (e.error()) {
    case EINVAL:
      // TODO: check if the server has been stopped
      break;
    default:
      printError(e.what());
    }
  } catch (const std::exception & e) {
    printError(e.what());
  } catch (...) {
    cerr << "Unknown connection error." << '\n';
    syslog(LOG_ERR, "Unknown connection error.");
  }
}

//----------------------------------
// ThreadedServer::ConnectionHandler
//----------------------------------

ThreadedServer::ConnectionHandler::ConnectionHandler(
    ThreadedServer         & server,
    shared_ptr<ITransport>   transport,
    shared_ptr<IProtocol>    protocol,
    IProcessor             & processor) :
    // initializers
    m_server    (server),
    m_transport (transport),
    m_protocol  (protocol),
    m_processor (processor) {
}

void ThreadedServer::ConnectionHandler::operator () () {
  try {
    for (;;) // until an exception is thrown
      m_processor.process(*m_protocol);
  } catch (const EofException &) {
    // connection closed, no need for alarm
  } catch (const std::exception & e) {
    printError(e.what());
  } catch (...) {
    cerr << "Unknown error." << '\n';
    syslog(LOG_ERR, "Unknown error.");
  }
  m_protocol.reset();
  m_transport.reset();
}
