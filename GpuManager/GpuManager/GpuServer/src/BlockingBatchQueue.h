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
class BlockingBatchQueue : public IBlockingQueue<T>
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
  BlockingBatchQueue() : m_interrupted(false) {}

  virtual void push(const T & item) {
    {
      scoped_lock lock(m_mutex);
      m_queue.push_back(item);
    }
    m_condition.notify_one();
  }

  virtual void pop(std::string & name, std::vector<T> & batch) {
    scoped_lock lock(m_mutex);

    while (m_queue.empty() && !m_interrupted)
      m_condition.wait(lock);
    if (m_interrupted)
      throw typename IBlockingQueue<T>::InterruptedError();

    name = m_queue.front()->Name();

    // move all items of this name from the queue and into the batch
    typename Queue::iterator p = stable_partition(m_queue.begin(), m_queue.end(), NotHasName(name));
    batch.assign(p, m_queue.end());
    m_queue.erase(p, m_queue.end());
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

  boost::condition_variable m_condition;
  boost::mutex              m_mutex;
};
