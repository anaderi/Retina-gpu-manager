#include "DataSender.h"

#include "GpuIpc/LocalSocketClient.h"
#include "GpuIpc/Protocol.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdint.h>
#include <thread>

using namespace boost;
using namespace boost::filesystem;
using namespace std;
using namespace std::chrono;

namespace {
  void readStream(ifstream & stream, void * data, size_t size) {
    stream.read(reinterpret_cast<char *>(data), size);
  }
}

void DataSender::print(const char * msg) const {
  ostringstream str;
  str << m_index << "-" << m_iteration << ": "<< msg << "\n";
  cout << str.str();
  cout.flush();
}

DataSender::DataSender(
    int                      index,
    const char             * servicePath,
    directory_entry_vector & paths,
    vector<DiffMessage>    & diffMessages,
    mutex                  & pathsMutex,
    bool                     verifyOutput,
    PerfLog                & perfLog) :
    m_verifyOutput (verifyOutput),
    m_servicePath  (servicePath),
    m_index        (index),
    m_paths        (paths),
    m_diffMessages (diffMessages),
    m_perfLog      (perfLog),
    m_mutex        (pathsMutex) {
}

void DataSender::operator() ()
try {
  LocalSocketClient transport (m_servicePath.c_str());
  Protocol          protocol  (transport);

  for (m_iteration = 0;; ++m_iteration) {
    // load and send the next item from the paths collection
    string path;
    {
      scoped_lock lock(m_mutex);
      if (m_paths.empty())
        break;
      path = m_paths.back().path().string();
      m_paths.pop_back();
    }

    auto s = steady_clock::now();

    string handlerName;
    vector<uint8_t> recordedInput;
    vector<uint8_t> recordedOutput;
    readData(path.c_str(), handlerName, recordedInput, recordedOutput);

    // send the name of the addressee
    protocol.writeString(handlerName);

    // send recorded input
    protocol.writeUInt32(recordedInput.size());
    if (!recordedInput.empty())
      protocol.writeData(recordedInput.data(), recordedInput.size());

    // receive the result and handle errors
    const size_t resultSize = protocol.readUInt32();

    const size_t FAIL_FLAG = 0xFFFFFFFF;
    if (resultSize == FAIL_FLAG) {
      string message = protocol.readString();
      ostringstream msg;
      msg << path << " | " << handlerName << " | error: " << message << "\n";
      cerr << msg.str();
      return;
    }

    vector<uint8_t> result(resultSize);
    if (!result.empty())
      protocol.readData(result.data(), resultSize);

    auto f = steady_clock::now();

    // record performance log
    {
      scoped_lock lock(m_mutex);
      m_perfLog.addRecord(time(0), duration_cast<duration<double>>(f - s).count());
    }

    // verify output
    if (m_verifyOutput) {
      string message = diff(result, recordedOutput);
      if (!message.empty()) {
        scoped_lock lock(m_mutex);
        m_diffMessages.push_back(DiffMessage(path, message));
      }
    }

    // receive performance information ­ seconds elapsed
    protocol.readDouble();
  }
} catch (const std::exception & e) {
  cout << "Unrecoverable error: " << e.what() << endl;
} catch (...) {
  cout << "Unkonwn unrecoverable error." << endl;
}

string DataSender::diff(
    const vector<uint8_t> & data,
    const vector<uint8_t> & reference) {
  // check size
  if (data.size() != reference.size()) {
    ostringstream message;
    message << "size mismatch: " << data.size() << " B instead of " << reference.size() << " B";
    return message.str();
  }

  // check contents
  size_t mismatchCount = 0;
  for (size_t i = 0, size = data.size(); i != size; ++i)
    mismatchCount += (data[i] != reference[i]);
  if (mismatchCount > 0) {
    ostringstream message;
    message << "data mismatch on " << mismatchCount << " bytes";
    return message.str();
  }

  return "";
}

void DataSender::readData(
    const char      * path,
    string          & handlerName,
    vector<uint8_t> & input,
    vector<uint8_t> & output) {
  ifstream stream(path, ios_base::binary);
  if (!stream) {
    ostringstream message;
    message << "File '" << path << "' could not be opened.";
    throw runtime_error(message.str());
  }

  uint32_t handlerSize;
  readStream(stream, &handlerSize, 4);

  // read name
  vector<char> handlerChars(handlerSize + 1); // +1 for terminating zero
  readStream(stream, handlerChars.data(), handlerSize);
  handlerName = handlerChars.data();

  // read input
  uint32_t inputSize;
  readStream(stream, &inputSize, 4);

  input.resize(inputSize);
  readStream(stream, input.data(), inputSize);

  // read output
  uint32_t outputSize;
  readStream(stream, &outputSize, 4);

  output.resize(outputSize);
  readStream(stream, output.data(), outputSize);
}
