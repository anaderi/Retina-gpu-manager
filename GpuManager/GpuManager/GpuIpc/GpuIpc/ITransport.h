#pragma once

#include <cstddef>

/// Transport layer interface.
class ITransport {
  public:
    virtual ~ITransport() {}

    /// Read a block of bytes.
    virtual void readBytes(void * data, size_t size) = 0;

    /// Write a block of bytes.
    virtual void writeBytes(const void * data, size_t size) = 0;
};
