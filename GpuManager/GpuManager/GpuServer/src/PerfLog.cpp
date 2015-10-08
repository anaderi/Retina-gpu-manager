#include "PerfLog.h"

using namespace std;

PerfLog::PerfLog(const char * filePath) :
    m_filePath (filePath) {
}

void PerfLog::addRecord(
    time_t       timestamp,
    const char * kernelName,
    double       seconds,
    size_t       inputBytes,
    size_t       outputBytes,
    size_t       batchCount) {
  prepareStream();

  m_out << timestamp          << '\t';
  m_out << escape(kernelName) << '\t';
  m_out << seconds            << '\t';
  m_out << inputBytes         << '\t';
  m_out << outputBytes        << '\t';
  m_out << batchCount         << '\n';
  m_out << flush;
}

string PerfLog::escape(const char * str) {
  string result(str);
  for (size_t i = 0, size = result.size(); i != size; ++i) {
    switch (str[i])
    {
      case '\t': case '\r': case '\n':
        result[i] = ' ';
    }
  }
  return result;
}

void PerfLog::prepareStream() {
  // lazy log file use
  if (!m_out.is_open()) {
    m_out.open(m_filePath.c_str(), ofstream::out | ofstream::app);

    // append header, if file is empty
    if (m_out.tellp() == 0)
      m_out << "timestamp\tname\ttime(s)\tin-size(B)\tout-size(B)\tcount\n" << flush;
  }
}
