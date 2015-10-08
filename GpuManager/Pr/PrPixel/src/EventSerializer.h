#pragma once

#include "PrPixelHit.h"
#include "PrPixelTrack.h"

#include <ostream>
#include <vector>

class EventSerializer {
public:
  typedef std::vector<uint8_t> Data;
  typedef std::vector<PrPixelHit*> HitMap;

public:
  EventSerializer();

  Data serializeEvent(PrPixelHit * hits, size_t hitCount);

  std::vector<PrPixelTrack> deserializeTracks(const Data & data) const;

  void printEvent(const Data & eventData, std::ostream & stream);

private:
  static size_t countModules(const PrPixelHit * hits, size_t hitCount);

  void reset();

private:
  HitMap m_hits;
  bool   m_hasData;
};
