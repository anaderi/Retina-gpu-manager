#pragma once

#include <vector>

#include <thread>
#include <memory>
#include <mutex>


/// Similar to boost::thread_group, but waits for thread completion on
/// destruction and automatically removes completed threads.
class ThreadGroup {
private:
  using lock_guard = std::lock_guard<std::mutex>;

public:
  ThreadGroup();
  ~ThreadGroup();

  /// Create a new std::thread object as-if by new thread(threadfunc) and add it
  /// to the group. 
  template<typename F>
  std::thread * createThread(F threadFunc) {
    using namespace std;

    lock_guard guard(m_mutex);

    // Call threadFunc, and then remove the thread from m_threads by thread ID.
    auto selfDeletingFunc = [=]() mutable {
      threadFunc();

      auto sameID = [](const thread * t) {
        return t->get_id() == this_thread::get_id();
      };

      lock_guard guard(m_mutex);
      thread::id id(this_thread::get_id());
      m_threads.erase(find_if(begin(m_threads), end(m_threads), sameID));
    };

    auto_ptr<thread> newThread(new thread(selfDeletingFunc));
    m_threads.push_back(newThread.get());
    return newThread.release();
  }

  /// Call join() on every running thread in the group.
  void joinAll();

  /// Returns the number of running threads in the group.
  std::size_t size();

private:
  std::vector<std::thread*> m_threads;
  std::mutex                m_mutex;
};
