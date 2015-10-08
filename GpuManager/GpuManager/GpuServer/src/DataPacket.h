#pragma once

#include "Data.h"

#include <string>

#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>

class DataPacket {
  private:

    typedef boost::mutex::scoped_lock scoped_lock;

  public:

    DataPacket(const std::string & name, const Data * data, Data * result);

    const std::string & Name()    { return m_name;    }
    const Data      * Input()   { return m_input;   }
    Data            * Output()  { return m_output;  }

    bool        ExceptionThrown()  { return m_exceptionThrown;  }
    std::string ExceptionMessage() { return m_exceptionMessage; }

    void SetExceptionMessage(const std::string & message);

    void Signal();
    void Wait();

  private:

    std::string    m_name;
    const Data * m_input;
    Data       * m_output;

    bool        m_exceptionThrown;
    std::string m_exceptionMessage;

    bool m_signalled;

    boost::condition_variable m_condition;
    boost::mutex              m_mutex;
};
