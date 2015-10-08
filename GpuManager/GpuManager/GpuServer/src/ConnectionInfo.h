#pragma once

#include <string>

// Simulating algebraic types for ConnectionInfo.
// The intention is to have a well-defined number of subtypes, enumerated in
// ConnectionInfo::Type. This makes it safe for the user to upcast from
// ConnectionInfo.

class ConnectionInfo {
	public:
		enum Type { Local, Tcp };

	protected:
		ConnectionInfo(Type type) : m_type(type) {}
		virtual ~ConnectionInfo() {}

	public:
		Type type() const { return m_type; }

	private:
		Type m_type;
};

class LocalConnectionInfo : public ConnectionInfo {
	public:
		LocalConnectionInfo(const char * path) : ConnectionInfo(Local), m_path(path) {}

		const char * path() const { return m_path.c_str(); }

		std::string m_path;
};

class TcpConnectionInfo : public ConnectionInfo {
	public:
		TcpConnectionInfo(const char * host, int port) : ConnectionInfo(Tcp), m_host(host), m_port(port) {}

		const char * host() const { return m_host.c_str(); }
		int          port() const { return m_port;         }

	private:
		std::string m_host;
		int         m_port;
};
