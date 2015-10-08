#pragma once

#include "IProtocol.h"

/// Request processor. Must be thread-safe.
class IProcessor {
  public:
    /// Process a request.
    virtual void process(IProtocol & protocol) = 0;
};
