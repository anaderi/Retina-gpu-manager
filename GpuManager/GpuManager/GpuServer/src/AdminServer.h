#pragma once

#include <map>
#include <set>
#include <stdexcept>
#include <stdint.h>
#include <vector>
#include "GpuIpc/IProcessor.h"

class IApp;
class IProtocol;

class AdminServer : public IProcessor {
  public: // interface
    AdminServer(IApp & app);

  public:
    virtual void process(IProtocol & protocol);

    void processExit        (IProtocol & protocol);
    void processLoadHandler (IProtocol & protocol);

  private:
    void exit();

  private:

  private:
    IApp & m_app;
};
