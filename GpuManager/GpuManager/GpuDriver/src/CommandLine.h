#pragma once

#include <string>

class CommandLine {
  public:

    CommandLine(const char * servicePath);

    bool parse(int argc, char * argv[]);

    std::string dataPath()     const { return m_dataPath;     }
    std::string servicePath()  const { return m_servicePath;  }
    int         threadCount()  const { return m_threadCount;  }
    bool        verifyOutput() const { return m_verifyOutput; }

  private:

    std::string m_dataPath;
    std::string m_servicePath;
    int         m_threadCount;
    bool        m_verifyOutput;
};
