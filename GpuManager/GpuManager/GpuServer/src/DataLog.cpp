#include "DataLog.h"

#include <boost/filesystem.hpp>

#include <sstream>
#include <fstream>
#include <stdint.h>

using namespace std;
using namespace boost::filesystem;

namespace {
  void writeStream(ofstream & stream, const void * data, size_t size) {
    stream.write(reinterpret_cast<const char *>(data), size);
  }

  string makePath(const string & dir, int counter) {
    ostringstream path;
    path << dir << path::preferred_separator << counter << ".dat";
    return path.str();
  }
}

DataLog::DataLog(bool enabled, const string & dir) :
    m_enabled (enabled),
    m_dir     (dir),
    m_counter (0) {
}

void DataLog::addRecord(
    const string & handlerName,
    const Data   & input,
    const Data   & output) {
  if (!m_enabled)
    return;

  create_directory(m_dir);

  string   path   (makePath(m_dir, m_counter++));
  ofstream stream (path.c_str(), ios::binary | ios::trunc);
  if (!stream)
    return;

  const uint32_t handlerSize = handlerName.size();
  writeStream(stream, &handlerSize, 4);
  writeStream(stream, handlerName.c_str(), handlerSize);

  const uint32_t inputSize = input.size();
  writeStream(stream, &inputSize, 4);
  writeStream(stream, input.data(), inputSize);

  const uint32_t outputSize = output.size();
  writeStream(stream, &outputSize, 4);
  writeStream(stream, output.data(), outputSize);
}

