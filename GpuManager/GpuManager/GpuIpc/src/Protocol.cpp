#include "Protocol.h"

#include "ITransport.h"

#include <vector>

using namespace std;

//-----------------
// public interface
//-----------------

Protocol::Protocol(ITransport & transport) :
    m_transport(transport) {
}

Protocol::~Protocol() {
}

//-------------------------
// IProtocol implementation
//-------------------------

bool Protocol::readBool() {
  uint8_t b = false;
  m_transport.readBytes(reinterpret_cast<void*>(&b), 1);
  return b;
}

void Protocol::writeBool(bool b) {
  m_transport.writeBytes(reinterpret_cast<void*>(&b), 1);
}

double Protocol::readDouble() {
  double n = 0.0;
  m_transport.readBytes(reinterpret_cast<void*>(&n), 8);
  return n;
}

void Protocol::writeDouble(double n) {
  m_transport.writeBytes(reinterpret_cast<void*>(&n), 8);
}

uint32_t Protocol::readUInt32() {
  uint32_t n = 0;
  m_transport.readBytes(reinterpret_cast<void*>(&n), 4);
  return n;
}

void Protocol::writeUInt32(uint32_t n) {
  m_transport.writeBytes(reinterpret_cast<void*>(&n), 4);
}

string Protocol::readString() {
  uint32_t length = readUInt32();
  vector<char> text(length + 1); // + terminating zero
  m_transport.readBytes(&text[0], length);
  return &text[0];
}

void Protocol::writeString(const string & str) {
  writeUInt32(str.size());
  m_transport.writeBytes(str.c_str(), str.size());
}

void Protocol::readData(void * data, size_t size) {
  m_transport.readBytes(data, size);
}

void Protocol::writeData(const void * data, size_t size) {
  m_transport.writeBytes(data, size);
}
