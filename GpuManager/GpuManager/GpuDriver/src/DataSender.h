#pragma once

#include "PerfLog.h"

#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/thread/mutex.hpp>

class DataSender {
  public:

    struct DiffMessage {
      std::string path;
      std::string message;

      DiffMessage(const std::string & path, const std::string & message)
        : path(path), message(message)
      {}
    };

    typedef std::vector<boost::filesystem::directory_entry> directory_entry_vector;

  private:

    typedef boost::mutex::scoped_lock scoped_lock;

  public:

    DataSender(
        int                        index,
        const char               * servicePath,
        directory_entry_vector   & paths,
        std::vector<DiffMessage> & diffMessages,
        boost::mutex             & pathsMutex,
        bool                       verifyOutput,
        PerfLog                  & perfLog);

    void operator() ();

  private:

    std::string diff(
        const std::vector<uint8_t> & data,
        const std::vector<uint8_t> & reference);

    static void readData(
      const char           * path,
      std::string          & handlerName,
      std::vector<uint8_t> & input,
      std::vector<uint8_t> & output);

    // helper for debugging purposes
    void print(const char * msg) const;

  private:

    const bool        m_verifyOutput;
    const std::string m_servicePath;

    int                        m_index;
    directory_entry_vector   & m_paths;
    std::vector<DiffMessage> & m_diffMessages;
    PerfLog                  & m_perfLog;
    boost::mutex             & m_mutex;

    int m_iteration;
};
