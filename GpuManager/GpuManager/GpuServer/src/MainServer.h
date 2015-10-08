#pragma once

#include "IBlockingQueue.h"
#include "DataPacket.h"
#include "GpuIpc/IProcessor.h"

#include "GpuHandler/IGpuHandler.h"

#include <map>
#include <memory>
#include <stdint.h>
#include <string>
#include <vector>

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

class IProtocol;
class DataLog;
class PerfLog;

class MainServer : public IProcessor {
private:

  using HandlerMap  = std::map<std::string, IGpuHandler*>;
  using IQueue      = IBlockingQueue<DataPacket*>;
  using scoped_lock = boost::mutex::scoped_lock;

public: // interface

  MainServer(PerfLog & perfLog, DataLog & dataLog, size_t batchSize);
  ~MainServer();

  void start();
  void stop();

  void loadHandler(const std::string & handlerName);

public: // IProcess implementation

  virtual void process(IProtocol & protocol);

private: // private functions

  std::string createInvalidHandlerMsg(const std::string & handler) const;

  IGpuHandler * getHandlerByName(const std::string & name);

  void processQueue();

private:

  // about thread safety:
  //  - m_handlers is protected by m_mutex
  //  - m_perfLog is only accessed serially
  //  - m_dataLog is thread-safe
  //  - m_dataQueue is thread-safe

  HandlerMap m_handlers;

  PerfLog & m_perfLog;
  DataLog & m_dataLog;

  std::unique_ptr<IQueue> m_dataQueue;

  boost::thread m_processingThread;

  boost::mutex m_mutex;
};
