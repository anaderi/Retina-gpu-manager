#pragma once

#include <fstream>
#include <string>

class PerfLog {
  public:

    PerfLog(const char * filePath);

    void addRecord(time_t timestamp, double secods);

  private:

    void prepareStream();

  private:

    std::string   m_filePath;
    std::ofstream m_out;
};
