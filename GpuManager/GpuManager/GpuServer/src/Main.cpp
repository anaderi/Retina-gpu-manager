#include "App.h"
#include "DataLog.h"
#include "CommandLine.h"
#include "ConnectionInfo.h"
#include "Controller.h"
#include "PerfLog.h"
#include "Logger.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <syslog.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

//#include <Gaudi/PluginService.h>

using namespace boost;
using namespace std;

// Run the process as a daemon.
void doDaemonize() {
  // fork off the parent process
  pid_t pid(fork());
  if (pid > 0)
    exit(EXIT_SUCCESS);
  if (pid < 0)
    throw runtime_error("Fork failed.");

  // change the file mode mask
  umask(0);

  // create a new SID for the child process
  if (setsid() < 0)
    throw runtime_error("Failed to create a new child process SID.");

  // change the current working directory
  if (chdir("/") < 0)
    throw runtime_error("Failed to change the current working directory.");
}

std::shared_ptr<ConnectionInfo> makeConnectionInfo(
    const char * connection,
    const char * localPath,
    const char * host,
    int          port) {
  if (strcmp(connection, "local") == 0) return std::make_shared<LocalConnectionInfo>(localPath);
  if (strcmp(connection, "tcp")   == 0) return std::make_shared<TcpConnectionInfo>(host, port);
  throw runtime_error("Invalid connection type. Must be one of: local, tcp.");
}

void exitPreviousInstance(Logger & logger, const char * path)
try {
  Controller controller(logger, path);
  controller.stopServer();
} catch (const std::exception & e) {
  logger.printError(e.what());
  // it's ok
}

bool anotherInstanceExists(Logger & logger, const char * path)
try {
  Controller controller(logger, path);
  return true;
} catch (const std::exception &) {
  return false;
}

// Main entry point.
int main(int argc, char * argv[])
try {
  const char * const defaultPath = "/tmp/GpuManager";
  const char * const defaultHost = "localhost";
  const int          defaultPort = 65000;

  CommandLine cl(defaultPath, defaultHost, defaultPort);
  if (!cl.parse(argc, argv))
    return EXIT_SUCCESS;

  if (cl.daemonize())
    doDaemonize();

  const bool useStdIO = !cl.daemonize();
  Logger logger(useStdIO);

  PerfLog perfLog("perf.log");

  const bool recordData = cl.dataDir()[0] != '\0';
  DataLog dataLog(recordData, cl.dataDir());

  const string adminPath(string(cl.localPath()) + "-admin");

  if (cl.exit()) {
    exitPreviousInstance(logger, adminPath.c_str());
    return EXIT_SUCCESS;
  }

  if (cl.handlerToLoad()[0] != '\0') {
    Controller controller(logger, adminPath.c_str());
    controller.loadHandler(cl.handlerToLoad());
    return EXIT_SUCCESS;
  }

  if (anotherInstanceExists(logger, adminPath.c_str())) {
    logger.printMessage("Another instance is still running. Terminating.");
    return EXIT_SUCCESS;
  }

  auto connectionInfo = makeConnectionInfo(cl.connectionType(), cl.localPath(), cl.host(), cl.port());
  App app(logger, perfLog, dataLog, adminPath.c_str(), cl.batchSize(), *connectionInfo);
  app.run();

  return EXIT_SUCCESS;
} catch (const std::exception & e) {
  cout << "Unrecoverable error: " << e.what() << '\n';
  syslog(LOG_ERR, "Unrecoverable error: %s", e.what());
  return EXIT_FAILURE;
} catch (...) {
  cout << "Unknown unrecoverable error." << '\n';
  syslog(LOG_ERR, "Unknown unrecoverable error.");
  return EXIT_FAILURE;
}
