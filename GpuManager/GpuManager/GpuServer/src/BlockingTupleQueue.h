#pragma once

#include "IBlockingQueue.h"

#include <algorithm>
#include <deque>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/thread/condition_variable.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

// requires:
// void (*T)::Name()

template <typename T>
class BlockingTupleQueue : public IBlockingQueue<T>
{
private:

  typedef boost::mutex::scoped_lock scoped_lock;

  typedef std::deque<T> Queue;

  struct NotHasName {
    const std::string & name;
    NotHasName(const std::string & name) : name(name) {}
    bool operator() (const T & item) { return item->Name() != name; }
  };

public:

  BlockingTupleQueue(size_t nToPop)
    : m_interrupted(false), m_nToPop(nToPop) {}

  virtual void push(const T & item) {
    scoped_lock lock(m_mutex);

    m_queue.push_back(item);

    if (m_queue.size() % m_nToPop == 0) {
      lock.unlock();
      m_condition.notify_one();
    }
  }

  virtual void pop(std::string & name, std::vector<T> & batch) {
    scoped_lock lock(m_mutex);

    while (m_queue.size() < m_nToPop && !m_interrupted)
      m_condition.wait(lock);
    if (m_interrupted)
      throw typename IBlockingQueue<T>::InterruptedError();

    name = m_queue.front()->Name();

    // move all items of this name from the queue and into the batch
    batch.assign(m_queue.end() - m_nToPop, m_queue.end());
    m_queue.erase(m_queue.end() - m_nToPop, m_queue.end());
  }

  virtual void interrupt() {
    {
      scoped_lock lock(m_mutex);
      m_interrupted = true;
    }
    m_condition.notify_all();
  }

private:

  Queue m_queue;

  bool m_interrupted;

  const size_t m_nToPop;

  boost::condition_variable m_condition;
  boost::mutex              m_mutex;
};
