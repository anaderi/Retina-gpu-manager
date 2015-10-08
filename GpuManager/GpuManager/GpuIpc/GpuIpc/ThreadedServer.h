#pragma once

#include <memory>

#include <boost/function.hpp>
#include <boost/thread/thread.hpp>

class IConnector;
class IProcessor;
class IProtocol;
class ITransport;

/// Multithreading server that creates one thread per connection.
class ThreadedServer {
private: // helper types

  class ConnectionHandler {
  public:
    ConnectionHandler (
      ThreadedServer                & server,
      std::shared_ptr<ITransport>   transport,
      std::shared_ptr<IProtocol>    protocol,
      IProcessor                    & processor);
    void operator () ();
  private:
    ThreadedServer                & m_server;
    std::shared_ptr<ITransport>   m_transport;
    std::shared_ptr<IProtocol>    m_protocol;
    IProcessor                    & m_processor;
  };

public: // public interface

  /// Factory for instantiating protocols for multiple connections.
  typedef boost::function<std::shared_ptr<IProtocol> (ITransport&)> ProtocolFactory;

  /// Create a ThreadedServer instance.
  ThreadedServer(
    IConnector      & connector,
    ProtocolFactory   protocolFactory,
    IProcessor      & processor);

  /// Creates and starts server thread.
  boost::thread * serve();

  /// Halts the server thread.
  void stop();

private: // private functions

  void serveConnector();

  bool stopRequested();

private: // data

  int m_threadID;

  IConnector      & m_connector;
  ProtocolFactory   m_protocolFactory;
  IProcessor      & m_processor;
};
