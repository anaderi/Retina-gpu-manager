#pragma once

#include "FastTrackSvc/IFastTrackSvc.h"

#include <GaudiKernel/Service.h>

class ITransport;
class IProtocol;
class Tracker;

class FastTrackSvc :
    public         Service,
    public virtual IFastTrackSvc {
  public:
    FastTrackSvc(const std::string & name, ISvcLocator * sl);
    virtual ~FastTrackSvc();

    // IFastTrackSvc implementation
    virtual void searchByPair(const std::vector<Track> & data, std::vector<int8_t> & result);

    // Service implementation
    virtual StatusCode queryInterface(const InterfaceID & riid, void ** ppvIF);
    virtual StatusCode initialize();
    virtual StatusCode finalize();

  private:
    void cleanup();
    void initIO();

  private:
    ITransport * m_transport;
    IProtocol  * m_protocol;
    Tracker    * m_tracker;
};
