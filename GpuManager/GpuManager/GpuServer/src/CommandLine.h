#pragma once

#include <string>

class CommandLine {
  public:
    CommandLine(const char * localPath, const char * tcpHost, int tcpPort);

    bool parse(int argc, char * argv[]);

    bool daemonize() const { return m_daemonize; }
    bool exit()      const { return m_exit;      }

    const char * dataDir()        const { return m_dataDir.c_str();       }
    const char * handlerToLoad()  const { return m_handlerToLoad.c_str(); }
    const char * connectionType() const { return m_connection.c_str();    }
    const char * localPath()      const { return m_localPath.c_str();     }
    const char * host()           const { return m_host.c_str();          }
    int          port()           const { return m_port;                  }
    size_t       batchSize()      const { return m_batchSize;        }

  private:
    bool        m_daemonize;
    bool        m_exit;
    std::string m_dataDir;
    std::string m_handlerToLoad;
    std::string m_connection;
    std::string m_localPath;
    std::string m_host;
    int         m_port;
    size_t      m_batchSize;
};
