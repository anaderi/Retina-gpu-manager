#include "CommandLine.h"

#include <iostream>

#include <boost/program_options.hpp>

using namespace boost::program_options;
using namespace std;

CommandLine::CommandLine(const char * servicePath) :
    m_servicePath  (servicePath),
    m_verifyOutput (false) {
}

bool CommandLine::parse(int argc, char * argv[]) {
  options_description desc("Supported options");
  desc.add_options()
    ("help", "display this help message"                                                              )
    ("service", value(&m_servicePath)->default_value(m_servicePath), "socket path"                    )
    ("threads", value(&m_threadCount)->default_value(1), "the number of requests to run concurrently" )
    ("verify",  value(&m_verifyOutput)->zero_tokens(),   "verify received output"                     )
    ("data",    value(&m_dataPath)->required(),      "a single gpuserver input record or a directory" );

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

  if (m_threadCount <= 0) {
    cout << desc << '\n';
    return false;
  }

  return true;
}
