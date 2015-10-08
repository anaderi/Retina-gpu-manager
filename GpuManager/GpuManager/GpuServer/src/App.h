#pragma once

#include "IApp.h"

#include "Logger.h"

#include "AdminServer.h"
#include "ConnectionInfo.h"
#include "MainServer.h"

#include "GpuIpc/Protocol.h"
#include "GpuIpc/TcpSocketServerConnector.h"
#include "GpuIpc/LocalSocketServerConnector.h"
#include "GpuIpc/ThreadedServer.h"

#include <memory>

class DataLog;
class PerfLog;

/** @class App App.h GpuServer/App.h
 * Application-wide API.
 *
 * @author Alexey Badalov
 * @date   2013-01-11
 */
class App : public IApp {
  public: // interface

    /// @param useStdIO switches between the use of the standard streams and
    /// syslog.
    App(
      Logger     & logger,
      PerfLog    & perfLog,
      DataLog    & dataLog,
      const char * adminPath,
      size_t       batchSize,
      const ConnectionInfo & connection);

    void run();

  public: // IApp implementation

    /// Safely terminates the application.
    virtual void exit();

    /// Load a GPU handler component.
    virtual void loadHandler(const std::string & handlerHanel);

  private: // private functions

    static std::shared_ptr<IProtocol>  getProtocol(ITransport & transport);
    static std::shared_ptr<IConnector> createConnector(const ConnectionInfo & info);

  private: // data

    Logger  & m_logger;

    LocalSocketServerConnector  m_adminConnector;
    std::shared_ptr<IConnector> m_mainConnector;

    ThreadedServer m_adminHost;
    ThreadedServer m_mainHost;

    MainServer  m_mainServer;
    AdminServer m_adminServer;
};
