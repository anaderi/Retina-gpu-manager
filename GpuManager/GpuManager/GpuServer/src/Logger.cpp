#include "Logger.h"

#include <cstdio>
#include <syslog.h>

Logger::Logger(bool useStdIO) :
    m_useStdIO (useStdIO) {
  if (!m_useStdIO) {
    freopen("/dev/null", "r", stdin);
    freopen("/dev/null", "w", stdout);
    freopen("/dev/null", "w", stderr);

    openlog("fasttrackd", LOG_CONS, LOG_DAEMON);
  }
}

Logger::~Logger() {
  if (m_useStdIO)
    closelog();
}

void Logger::printError(const char * text) {
  if (m_useStdIO)
    fprintf(stderr, "%s\n", text);
  else
    syslog(LOG_ERR, "%s\n", text);
}

void Logger::printMessage(const char * text) {
  if (m_useStdIO)
    fprintf(stdout, "%s\n", text);
  else
    syslog(LOG_INFO, "%s\n", text);
}
