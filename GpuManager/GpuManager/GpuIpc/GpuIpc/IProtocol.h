#pragma once

#include <stdint.h>
#include <string>
#include <vector>

/// Communications protocol layer interface.
class IProtocol {
  public:
    virtual ~IProtocol() {}

    /// Reads a boolean.
    virtual bool readBool() = 0;

    /// Writes a boolean.
    virtual void writeBool(bool b) = 0;

    /// Reads a double.
    virtual double readDouble() = 0;

    /// Writes a double.
    virtual void writeDouble(double n) = 0;

    /// Reads an unsigned 32-bit integer.
    virtual uint32_t readUInt32() = 0;

    /// Writes an unsigned 32-bit integer.
    virtual void writeUInt32(uint32_t n) = 0;

    /// Reads a standard string.
    virtual std::string readString() = 0;

    /// Writes a standard string.
    virtual void writeString(const std::string & str) = 0;

    /// Reads a block of data.
    virtual void readData(void * data, size_t size) = 0;

    /// Writes a block of data.
    virtual void writeData(const void * data, size_t size) = 0;
};
