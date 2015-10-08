#pragma once

class Logger {
  public: // interface

    Logger(bool useStdIO);
    ~Logger();

    /// Prints an error message.
    /// @param text is the message to be printed.
    void printError(const char * text);
    //
    /// Prints an information message.
    /// @param text is the message to be printed.
    void printMessage(const char * text);

  private: // data

    const bool m_useStdIO;
};
