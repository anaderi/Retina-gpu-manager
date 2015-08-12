#include "Tools.h"

#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>

EventInfo parseEvent(const uint8_t * input, size_t size)
{
  const uint8_t * end = input + size;
  int h_no_sensors       = *((int32_t*)input); input += sizeof(int32_t);
  int h_no_hits          = *((int32_t*)input); input += sizeof(int32_t);
  int* h_sensor_Zs        = ((int32_t*)input); input += sizeof(int32_t) * h_no_sensors;
  int* h_sensor_hitStarts = ((int32_t*)input); input += sizeof(int32_t) * h_no_sensors;
  int* h_sensor_hitNums   = ((int32_t*)input); input += sizeof(int32_t) * h_no_sensors;
  uint32_t* h_hit_IDs          = (uint32_t*)input; input += sizeof(uint32_t) * h_no_hits;
  float * h_hit_Xs           = (float*)  input; input += sizeof(float)   * h_no_hits;
  float * h_hit_Ys           = (float*)  input; input += sizeof(float)   * h_no_hits;
  float * h_hit_Zs           = (float*)  input; input += sizeof(float)   * h_no_hits;

  if (input != end)
    throw std::runtime_error("failed to deserialize event"); 

  EventInfo event;
  
  event.minX = *std::min_element(h_hit_Xs, h_hit_Xs + h_no_hits);
  event.minY = *std::min_element(h_hit_Ys, h_hit_Ys + h_no_hits);
  event.minZ = *std::min_element(h_hit_Zs, h_hit_Zs + h_no_hits);
  event.maxX = *std::max_element(h_hit_Xs, h_hit_Xs + h_no_hits);
  event.maxY = *std::max_element(h_hit_Ys, h_hit_Ys + h_no_hits);
  event.maxZ = *std::max_element(h_hit_Zs, h_hit_Zs + h_no_hits);
  
  {
    std::vector<Hit> parsedHits;
    parsedHits.reserve(h_no_hits);
    int count = 0;
    for (uint32_t sensorId =  0; sensorId < h_no_sensors; ++sensorId)
    {
      for (int i = 0; i < h_sensor_hitNums[sensorId]; ++i) 
      {
        int hitId = h_sensor_hitStarts[sensorId] + i;
        assert(hitId >= 0 && hitId < h_no_hits );
        parsedHits.emplace_back(
          h_hit_Xs[hitId],
          h_hit_Ys[hitId],
          h_hit_Zs[hitId],
          //h_hit_IDs[hitId],
          count++,
          sensorId
        );
      }
    }
    event.hits = parsedHits;
  }
  return event;
}

std::vector<uint8_t> putTracksInOutputFormat(
  const std::vector<Hit>& hits,
  const std::vector<Track>& tracks
)
{
  std::map<int, int> remap;
  for (int i = 0; i < hits.size(); i++)
  {
    remap[hits[i].id] = i;
  }
  std::vector<uint8_t> output(tracks.size() * sizeof(Track));
  Track * outputPtr = (Track*)&output[0];
  for (int j = 0; j < tracks.size(); j++)
  {
    /*for (int i = 0; i < tracks[j].hitsNum; i++)
    {
      auto hit = hits[tracks[j].hits[i]];
      std::cerr << "x: " << ((hit.x) * (_maxX - _minX) + _minX)
        << ", y: " << ((hit.y) * (_maxY - _minY) + _minY)
        << "z: " << ((hit.z) * (_maxZ - _minZ) + _minZ) << std::endl;
    }
     */
    for (int i = 0; i < tracks[j].hitsNum; i++)
    {
      auto a = hits[tracks[j].hits[i]];
      //fprintf(stderr, "[%f,%f,%f,%d],\n", a.x, a.y, a.z, j);      
    }
    //std::cerr << std::endl;
     
  }
  for (size_t i = 0; i != tracks.size(); ++i) 
  {
    Track copy = tracks[i];
    /*for (int j = 0; j < copy.hitsNum; j++)
    {
      copy.hits[j] = remap[copy.hits[j]];
    }*/
    outputPtr[i] = copy;
  }
  return output;
}

/*void printHit(const Hit& hit, std::ostream& stream) 
{
  stream << " " << std::setw(8) << hit.id 
    << " module " << std::setw(2) << hit.sensorId
    << ", x " << std::setw(6) << hit.x
    << ", y " << std::setw(6) << hit.y
    << ", z " << std::setw(6) << hit.z 
    << std::endl;
  
}

void printSolution(
  const std::vector<Track>& tracks, 
  const std::vector<Hit>& hits,
  std::ostream& stream
)
{
  std::map<uint32_t, Hit> hitMap;
  for (const Hit& hit: hits) {
    hitMap[hit.id] = hit;
    assert(hitMap[hit.id].sensorId == hit.sensorId);
  }
  int id = 0;
  for (size_t trackId = 0; trackId < tracks.size(); ++trackId) 
  {
    const Track& track = tracks[trackId];
    if (track.hitsNum != 0)
    {
      stream << "Track #" << id++ << ", length " << track.hitsNum << std::endl;
      for (int i = track.hitsNum - 1; i >= 0; --i)
      {
        printHit(hitMap[track.hits[i]], stream);
      }
      stream << std::endl;
    }
  }

  stream << std::endl;
}
*/
double getQuatile(std::vector<double> data, double ratio)
{
  std::sort(data.begin(), data.end());
  return data[data.size() * ratio];
}

