#pragma once

#include "IProtocol.h"

class ITransport;

/// The standard protocol implementation.
class Protocol : public IProtocol {
  public: // public interface

    /// Construct a protocol instance over a given transport layer.
    Protocol(ITransport & transport);

    virtual ~Protocol();

  public: // IProtocol implementation

    /// Reads a boolean.
    virtual bool readBool();

    /// Writes a boolean.
    virtual void writeBool(bool b);

    /// Reads a double.
    virtual double readDouble();

    /// Writes a double.
    virtual void writeDouble(double n);

    /// Reads an unsigned 32-bit integer.
    virtual uint32_t readUInt32();

    /// Writes an unsigned 32-bit integer.
    virtual void writeUInt32(uint32_t n);

    /// Reads a standard string.
    virtual std::string readString();

    /// Writes a standard string.
    virtual void writeString(const std::string & str);

    /// Reads a block of data.
    virtual void readData(void * data, size_t size);

    /// Writes a block of data.
    virtual void writeData(const void * data, size_t size);

  private:
    ITransport & m_transport;
};
