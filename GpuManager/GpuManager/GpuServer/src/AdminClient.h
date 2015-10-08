#pragma once

#include <string>

class IProtocol;

class AdminClient {
  public:
    AdminClient(IProtocol & protocol);

  public:
    void exit();
    void loadHandler(const std::string & handlerName);

  private:
    IProtocol & m_protocol;
};
