#include "CommandLine.h"

#include <iostream>
#include <stdexcept>

#include <boost/program_options.hpp>

using namespace boost::program_options;
using namespace std;

CommandLine::CommandLine(const char * localPath, const char * tcpHost, int tcpPort) :
    m_daemonize  (false),
    m_exit       (false),
		m_connection ("local"),
    m_localPath  (localPath),
		m_host       (tcpHost),
		m_port       (tcpPort) {
}

bool CommandLine::parse(int argc, char * argv[]) {
  options_description desc("Supported options");
  desc.add_options()
    ("help", "display this help message")
    ("daemonize",  bool_switch(&m_daemonize),                         "run the process as a daemon"     )
    ("exit",       bool_switch(&m_exit),                              "stop the server"                 )
    ("load",       value(&m_handlerToLoad),                           "load a handler component"        )
    ("connection", value(&m_connection)->default_value(m_connection), "one of: local, tcp"              )
    ("localPath",  value(&m_localPath)->default_value(m_localPath),   "local socket path"               )
    ("tcpHost",    value(&m_host)->default_value(m_host),             "tcp hostname or server IP"       )
    ("tcpPort",    value(&m_port)->default_value(m_port),             "tcp server port"                 )
    ("batchSize",  value(&m_batchSize)->default_value(0),             "positive value fixes batch size" )
    ("datadir",    value(&m_dataDir),                                 "save transactions to directory"  );

  variables_map vm;
  try {
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << '\n';
      return false;
    }
    notify(vm);
  } catch (const std::exception & e) {
    cerr << e.what() << endl;
    cout << desc << '\n';
    return false;
  }

  return true;
}
