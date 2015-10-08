#pragma once

#include <string>

class IApp {
  public:
    /// Safely terminates the application.
    virtual void exit() = 0;

    /// Load a GPU handler component.
    virtual void loadHandler(const std::string & handlerName) = 0;
};
