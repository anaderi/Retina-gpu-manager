#include "ThreadGroup.h"

using namespace std;

ThreadGroup::ThreadGroup() {
}

ThreadGroup::~ThreadGroup() {
  lock_guard guard(m_mutex);
  for (thread * t : m_threads) {
    if (t->joinable())
      t->join();
    delete t;
  }
}

void ThreadGroup::joinAll() {
  lock_guard guard(m_mutex);
  for (thread * t : m_threads) {
    if (t->joinable())
      t->join();
  }
}

size_t ThreadGroup::size() {
  lock_guard guard(m_mutex);
  return m_threads.size();
}
