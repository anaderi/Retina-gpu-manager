#include "PrPixelCudaHandler.h"
#include "Retina.h"

#include <algorithm>
#include <vector>

using namespace std;

DECLARE_COMPONENT(PrPixelCudaHandler)

void PrPixelCudaHandler::operator() (
    const Batch & batch,
    Alloc         allocResult,
    AllocParam    allocResultParam) {
  // gpuPixelSearchByTriplet handles several events in parallel
  vector<Data> trackCollection(batch.size());
  //gpuPixelSearchByTriplet(batch, trackCollection);
  cpuRetinaInvocation(batch, trackCollection);
 
  for (int i = 0, size = batch.size(); i != size; ++i) {
    const Data & t = trackCollection.at(i);
    copy(t.begin(), t.end(), allocResult(i, t.size(), allocResultParam));
  }
}
