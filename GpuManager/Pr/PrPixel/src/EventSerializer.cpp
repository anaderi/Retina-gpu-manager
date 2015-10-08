#include "EventSerializer.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>

#include <boost/functional/hash.hpp>

using namespace std;

using Data = EventSerializer::Data;

//----------------------------
// Serialized track structure.
//----------------------------

const size_t MAX_TRACK_SIZE = 24;
struct GpuTrack {
  uint32_t hitsNum;
  uint32_t hits[MAX_TRACK_SIZE];
};

//----------------------------
// EventSerializer definition.
//----------------------------

EventSerializer::EventSerializer() : m_hasData(false) {
}

Data EventSerializer::serializeEvent(PrPixelHit * hits, size_t hitCount) {
  if (m_hasData)
    reset();

  const uint32_t sensorCount = static_cast<uint32_t>(countModules(hits, hitCount));

  const size_t sensorCountSize     = sizeof(int32_t);
  const size_t hitCountSize        = sizeof(int32_t);
  const size_t sensorZsSize        = sensorCount * sizeof(int32_t);
  const size_t sensorHitStartsSize = sensorCount * sizeof(int32_t);
  const size_t sensorHitsNumsSize  = sensorCount * sizeof(int32_t);
  const size_t hitIDsSize          = hitCount    * sizeof(int32_t);
  const size_t hitXsSize           = hitCount    * sizeof(float);
  const size_t hitYsSize           = hitCount    * sizeof(float);
  const size_t hitZsSize           = hitCount    * sizeof(float);

  Data buffer(sensorCountSize + hitCountSize + sensorZsSize + sensorHitStartsSize
      + sensorHitsNumsSize + hitIDsSize + hitXsSize + hitYsSize + hitZsSize);

  uint8_t * dst = buffer.data();

  *(int32_t *)dst = sensorCount; dst += sensorCountSize;
  *(int32_t *)dst = hitCount;    dst += hitCountSize;

  int32_t * sensorZs          = (int32_t*) dst; dst += sensorZsSize;
  int32_t * sensorHitStarts   = (int32_t*) dst; dst += sensorHitStartsSize;
  int32_t * sensorHitNums     = (int32_t*) dst; dst += sensorHitsNumsSize;
  int32_t * hitIDs            = (int32_t*) dst; dst += hitIDsSize;
  float   * hitXs             = (float*)   dst; dst += hitXsSize;
  float   * hitYs             = (float*)   dst; dst += hitYsSize;
  float   * hitZs             = (float*)   dst; dst += hitZsSize;

  int lastModule = -1;

  if (dst != buffer.data() + buffer.size())
    throw runtime_error("event serialization failed");

  size_t sensorCountCheck = 0;

  for (size_t i = 0; i != hitCount; ++i) {
    PrPixelHit & hit = hits[i];

    if (lastModule != hit.module()) {
      lastModule = hit.module();
      *sensorHitStarts = i;       ++sensorHitStarts;
      *sensorHitNums   = 1;       ++sensorHitNums;
      *sensorZs        = hit.z(); ++sensorZs;
      ++sensorCountCheck;
    } else {
      ++*(sensorHitNums - 1);
    }

    *hitIDs = hit.id().vpID(); ++hitIDs;
    *hitXs  = hit.x();         ++hitXs;
    *hitYs  = hit.y();         ++hitYs;
    *hitZs  = hit.z();         ++hitZs;

    m_hits.push_back(&hit);
  }

  if (sensorCount != sensorCountCheck)
    throw runtime_error("sensor counting failed during event serialization");

  m_hasData = true;

  //printEvent(buffer, cout);

  return buffer;
}

vector<PrPixelHit*> mapHits(
    const GpuTrack                & track,
    const EventSerializer::HitMap & hitMap) {
  vector<PrPixelHit*> hits;
  for (size_t i = 0; i != track.hitsNum; ++i)
    hits.push_back(hitMap.at(track.hits[i]));
  return hits;
}

vector<PrPixelTrack> EventSerializer::deserializeTracks(const Data & data) const {
  if (!m_hasData)
    throw runtime_error("deserialization using an empty serializer");

  vector<PrPixelTrack> tracks;

  if (data.empty())
    throw runtime_error("empty track data");
  if (data.size() % sizeof(GpuTrack) != 0) {
    ostringstream msg;
    msg << "invalid track data size: " << data.size();
    throw runtime_error(msg.str());
  }

  const GpuTrack * gpuTracks = reinterpret_cast<const GpuTrack*>(data.data());

  const size_t count = data.size() / sizeof(GpuTrack);

  //print(gpuTracks, count);

  for (size_t i = 0; i != count; ++i) {
    const GpuTrack & t = gpuTracks[i];

    if (t.hitsNum < 1 || t.hitsNum > MAX_TRACK_SIZE) {
      ostringstream msg;
      ostringstream() << "invalid track hitsNum: " << t.hitsNum;
      throw runtime_error(msg.str());
    }
    if (i == 0)
    {
      for (size_t i = 0; i != t.hitsNum; ++i) {
        auto hit = *m_hits.at(t.hits[i]);
        std::cerr << "x: " << hit.x() << ", y: " << hit.y() << "z: " << hit.z() << std::endl;
      }
    }
    tracks.push_back(PrPixelTrack(mapHits(t, m_hits)));
  }
  //std::cerr << tracks[0].m_hits[0]->x() << std::endl;
  return tracks;
}

template <typename T>
size_t hash(const T * p, size_t n) {
  return boost::hash_value(vector<T>(p, p + n));
}

void EventSerializer::printEvent(const Data & eventData, ostream & stream) {
  const uint8_t * data = eventData.data();

  const int sensorCount = *(const int32_t*)data; data += sizeof(int32_t);
  const int hitCount    = *(const int32_t*)data; data += sizeof(int32_t);

  const int32_t * sensorZs        = (const int32_t*) data; data += sensorCount * sizeof(int32_t);
  const int32_t * sensorHitStarts = (const int32_t*) data; data += sensorCount * sizeof(int32_t);
  const int32_t * sensorHitsNums  = (const int32_t*) data; data += sensorCount * sizeof(int32_t);
  const int32_t * hitIDs          = (const int32_t*) data; data += hitCount    * sizeof(int32_t);
  const float   * hitXs           = (const float*)   data; data += hitCount    * sizeof(float);
  const float   * hitYs           = (const float*)   data; data += hitCount    * sizeof(float);
  const float   * hitZs           = (const float*)   data; data += hitCount    * sizeof(float);

  if (data != eventData.data() + eventData.size())
    throw runtime_error("invalid event format");

  stream << sensorCount << " sensors, " << hitCount << " hits" << endl;
  stream << "sensor zs: 0x" << hex << ::hash(sensorZs,        sensorCount) << dec << endl;
  stream << "sensor hs: 0x" << hex << ::hash(sensorHitStarts, sensorCount) << dec << endl;
  stream << "sensor hn: 0x" << hex << ::hash(sensorHitsNums,  sensorCount) << dec << endl;
  stream << "hit id:    0x" << hex << ::hash(hitIDs,          hitCount)    << dec << endl;
  stream << "hit xs:    0x" << hex << ::hash(hitXs,           hitCount)    << dec << endl;
  stream << "hit ys:    0x" << hex << ::hash(hitYs,           hitCount)    << dec << endl;
  stream << "hit zs:    0x" << hex << ::hash(hitZs,           hitCount)    << dec << endl;
}

size_t EventSerializer::countModules(const PrPixelHit * hits, size_t hitCount) {
  size_t sensorCount = 0;
  int    lastModule  = -1;
  for (size_t i = 0; i != hitCount; ++i) {
    if (lastModule != hits[i].module()) {
      lastModule = hits[i].module();
      ++sensorCount;
    }
  }
  return sensorCount;
}

void EventSerializer::reset() {
  m_hasData = false;
  m_hits.clear();
}
