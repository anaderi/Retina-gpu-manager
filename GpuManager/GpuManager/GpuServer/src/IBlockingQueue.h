#pragma once

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
class IBlockingQueue
{
public:
  struct InterruptedError : std::runtime_error {
    InterruptedError() : std::runtime_error("BlockingBatchQueue was interrupted.") { }
  };

public:

  virtual ~IBlockingQueue() { }

  virtual void push(const T & item) = 0;

  virtual void pop(std::string & name, std::vector<T> & batch) = 0;

  virtual void interrupt() = 0;
};
