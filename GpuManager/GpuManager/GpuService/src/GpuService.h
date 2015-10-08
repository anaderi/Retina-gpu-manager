#pragma once

#include "GpuService/IGpuService.h"

#include <GaudiKernel/Property.h>
#include <GaudiKernel/Service.h>

class ITransport;
class IProtocol;
class Tracker;

class GAUDI_API GpuService : public extends1<Service, IGpuService> {
public: // construction
  GpuService(const std::string & name, ISvcLocator * sl);
  virtual ~GpuService();

public: // IGpuService implementation
  virtual void submitData(
      std::string  handlerName,
      const void * data,
      const size_t size,
      Alloc        allocResults,
      AllocParam   allocResultsParam);

public: // Service implementation
  virtual StatusCode initialize();
  virtual StatusCode finalize();

private: // private functions
  void releaseIO();
  void acquireIO();

private:
  ITransport * m_transport;
  IProtocol  * m_protocol;

  StringProperty  m_connectionType;
  StringProperty  m_localSocketPath;
  StringProperty  m_serverHost;
  IntegerProperty m_serverPort;
};
