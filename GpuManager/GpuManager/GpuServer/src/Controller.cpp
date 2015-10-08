#include "Controller.h"

using namespace std;

Controller::Controller(
    Logger     & logger,
    const char * adminPath) :
    // initializers
    m_logger    (logger),
    m_transport (adminPath),
    m_protocol  (m_transport),
    m_admin     (m_protocol) {
}

void Controller::stopServer() {
  m_admin.exit();
}

void Controller::loadHandler(const string & handlerName) {
  m_admin.loadHandler(handlerName);
}
