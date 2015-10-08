#pragma once

#include <Gaudi/PluginService.h>

#include <stdint.h>
#include <string>
#include <vector>

// The handler is an object that takes a batch of input data sets
// and computes the results for each. This approach enables taking
// advantage of high-throughput high-latency processors.
class IGpuHandler
{
  public:

    virtual ~IGpuHandler();

    // Data is processed in batches. A batch consists of input data sets
    // from multiple clients.

    typedef std::vector<uint8_t>     Data;
    typedef std::vector<const Data*> Batch;

    // Allocator for the results of the compuation. It takes:
    // - the batch member index to which the result applies
    // - size of the memory block it returns
    // - an opaque parameter given to the handler

    typedef void * AllocParam;
    typedef uint8_t * (*Alloc)(
        std::size_t index,
        std::size_t size,
        AllocParam  param);

    // Main handler function. Using the allocator function allows the caller to
    // manage memory allocation and deallocation.

    virtual void operator() (
        const Batch & batch,
        Alloc         allocResult,
        AllocParam    allocResultParam) = 0;

    // Factory for plugin support.
    typedef Gaudi::PluginService::Factory<IGpuHandler*> Factory;
};
