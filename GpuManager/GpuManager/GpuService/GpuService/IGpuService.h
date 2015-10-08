#pragma once

#include "GaudiKernel/IInterface.h"

/// Base interface providing GPU-accelerated algorithms.
class IGpuService : public virtual IInterface {
  public: // helper types

    /// User-defined data to be passed through to the allocator.
    typedef void * AllocParam;

    /// Allocator function.
    typedef void * (*Alloc)(size_t size, AllocParam param);

  public: // interface

    virtual ~IGpuService() {}

    /// Submit a data package to a specific handler,
    /// supplying an allocator for the handler's result.
    virtual void submitData(
        std::string  handlerName,
        const void * data,
        const size_t size,
        Alloc        allocResults,
        AllocParam   allocResultsParam) = 0;

    /// Gaudi interface declaration
    DeclareInterfaceID(IGpuService, 1, 0);
};
