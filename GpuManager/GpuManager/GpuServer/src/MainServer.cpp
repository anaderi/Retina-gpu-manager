#include "BlockingBatchQueue.h"
#include "BlockingTupleQueue.h"
#include "DataLog.h"
#include "MainServer.h"
#include "PerfLog.h"

#include "GpuIpc/IProtocol.h"
#include "GpuHandler/IGpuHandler.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <ctime>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include <Gaudi/PluginService.h>

using namespace boost;
using namespace std;
using namespace std::chrono;

//--------
// helpers
//--------

namespace {
  template <typename T, size_t Size>
  constexpr size_t arraySize(const T (&)[Size]) {
    return Size;
  }

  IBlockingQueue<DataPacket*> * MakeDataQueue(size_t batchSize) {
    using T = DataPacket*;
    return batchSize > 0 ?
      static_cast<IBlockingQueue<T>*>(new BlockingTupleQueue<T>(batchSize)) :
      static_cast<IBlockingQueue<T>*>(new BlockingBatchQueue<T>());
  }

  uint8_t * allocVector(
      size_t index,
      size_t size,
      IGpuHandler::AllocParam param) {
    if (size == 0)
      return nullptr;
    using Batch = vector<Data*>;
    Batch * batch = reinterpret_cast<Batch*>(param);
    Data * data = batch->at(index);
    data->resize(size);
    return data->data();
  }
}

//----------
// interface
//----------

MainServer::MainServer(PerfLog & perfLog, DataLog & dataLog, size_t batchSize) :
    m_perfLog   (perfLog),
    m_dataLog   (dataLog),
    m_dataQueue (MakeDataQueue(batchSize)) {
}

MainServer::~MainServer()
{
  for (auto & i : m_handlers)
    delete i.second;
}

void MainServer::start() {
    m_processingThread = thread(&MainServer::processQueue, this);
}

void MainServer::stop() {
  m_dataQueue->interrupt();
}

void MainServer::loadHandler(const string & handlerName) {
  IGpuHandler * handler = IGpuHandler::Factory::create(handlerName);
  if (!handler) {
    ostringstream msg;
    msg << "could not load handler '" << handlerName << "'";
    throw runtime_error(msg.str());
  }
  {
    scoped_lock lock(m_mutex);

    // handle repeated loading of the same handler
    HandlerMap::const_iterator i = m_handlers.find(handlerName);
    if (i != m_handlers.end())
      delete i->second;

    m_handlers[handlerName] = handler;
  }
  cout << "loaded handler '" << handlerName << "'\n";
}

//------------------------
// IProcess implementation
//------------------------

void MainServer::process(IProtocol & protocol) {
  const size_t FAIL_FLAG = 0xFFFFFFFF;

  std::string handlerName = protocol.readString();

  size_t size = protocol.readUInt32();

  Data input(size);
  protocol.readData(&input[0], size);

  { // check whether the handler is loaded
    scoped_lock lock(m_mutex);
    if (m_handlers.find(handlerName) == m_handlers.end()) {
      // when the handler is not found, inform the client
      protocol.writeUInt32(FAIL_FLAG);
      protocol.writeString(createInvalidHandlerMsg(handlerName));
      return;
    }
  }

  // call the handler
  Data output;
  DataPacket packet(handlerName, &input, &output);
  m_dataQueue->push(&packet);

  auto s = steady_clock::now();
  packet.Wait();
  auto f = steady_clock::now();

  // forward exceptions to client
  if (packet.ExceptionThrown()) {
    protocol.writeUInt32(FAIL_FLAG);
    protocol.writeString(packet.ExceptionMessage());
    return;
  }

  m_dataLog.addRecord(handlerName, input, output);

  // return the output to client
  protocol.writeUInt32(output.size());
  protocol.writeData(&output[0], output.size());
  protocol.writeDouble(duration_cast<duration<double>>(f - s).count());
}

//------------------
// Private functions
//------------------

string MainServer::createInvalidHandlerMsg(const string & handler) const {
  ostringstream msg;
  msg << "invalid handler name: " << handler << "; ";
  msg << "valid handlers: ";
  bool first = true;
  for (const auto & i : m_handlers) {
    if (first)
      first = false;
    else
      msg << ", ";
    msg << i.first;
  }
  msg << ".";
  return msg.str();
}

IGpuHandler * MainServer::getHandlerByName(const std::string & name) {
  scoped_lock lock(m_mutex);
  HandlerMap::const_iterator i = m_handlers.find(name);
  if (i == m_handlers.end())
    throw runtime_error(createInvalidHandlerMsg(name));
  return i->second;
}

void MainServer::processQueue()
try {
  // the data queue throws an exception when interrupted
  while (true) {
    string              name;
    vector<DataPacket*> batch;
    m_dataQueue->pop(name, batch);

    IGpuHandler * handler = getHandlerByName(name);

    // prepare data
    vector<const Data*> input  (batch.size());
    vector<Data*>       output (batch.size());
    for (size_t i = 0, size = batch.size(); i != size; ++i) {
      input[i]  = batch[i]->Input();
      output[i] = batch[i]->Output();
    }

    double secondsElapsed = 0.0;
    try {
      // execute handler
      auto s = steady_clock::now();
      (*handler)(input, allocVector, &output);
      auto f = steady_clock::now();

      secondsElapsed = duration_cast<duration<double>>(f - s).count();
    } catch (const std::exception & e) {
      // propagate the exception to all client in the batch
      for (size_t i = 0, size = batch.size(); i != size; ++i)
        batch[i]->SetExceptionMessage(e.what());
    }

    // gather statistics
    auto addSize = [](size_t total, const Data * data) { return total + data->size(); };
    size_t totalInputSize  = accumulate(input.begin(),  input.end(),  0u, addSize);
    size_t totalOutputSize = accumulate(output.begin(), output.end(), 0u, addSize);
    m_perfLog.addRecord(
        time(0), name.c_str(), secondsElapsed,
        totalInputSize, totalOutputSize, batch.size());

    // wake up the clients
    // wait for them to finish before moving onto the next batch
    for (size_t i = 0, size = batch.size(); i != size; ++i)
      batch[i]->Signal();
  }
} catch (const IQueue::InterruptedError &) {
  // it's ok
  // someone wants us to terminate quickly
}
