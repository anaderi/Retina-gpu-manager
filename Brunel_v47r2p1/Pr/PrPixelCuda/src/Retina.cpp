#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <cstdlib>

#include "Retina.h"
#include "Grid.h"
#include "Physics.h"
#include "Tools.h"
#include "Logger.h"

TrackPure generateTrackFromIndex(
  const std::vector<int>& multiIndex,
  const std::vector<std::vector<double> > & dim
)
{
  return TrackPure(
    dim[0][multiIndex[0]],
    dim[1][multiIndex[1]],
    dim[2][multiIndex[2]],
    dim[3][multiIndex[3]]
  );
}

std::vector<TrackPure> restoreTracks(
  const std::vector<std::vector<double> >& dimensions,
  const std::vector<TrackPure>& grid,
  const std::vector<double>& responces
)
{
  int gridSizeNoBoarders = calculateGridSize(dimensions, 1);
  std::vector<int> indexes(dimensions.size(), 1);
  std::vector<TrackPure> restored;
  
  const double threshold = getQuatile(responces, TRACK_THRESHOLD);
  for (int i = 0; i < gridSizeNoBoarders; ++i)
  {
    std::vector<int> neighbours = generateNeighboursIndexes(indexes, dimensions);
    int currentIndex = multiIndexToIndex(indexes, dimensions);
    bool isLocalMaximum = true;
    for (int neighbour : neighbours)
    {
      if (responces[neighbour] > responces[currentIndex])
      {
        isLocalMaximum = false;
      }
    }
    isLocalMaximum = /*isLocalMaxumum ||*/ (responces[currentIndex] > threshold);
    if (isLocalMaximum)
    {
      TrackPure answer = grid[currentIndex] * responces[currentIndex];
      double sum_responce = responces[currentIndex];
      for (int neighbour : neighbours)
      {
        answer = answer + grid[neighbour] * responces[neighbour];
        sum_responce += responces[neighbour];
      }
      restored.push_back(answer * (1.0 / sum_responce));      
    }
    generateNextMultiIndex(indexes, dimensions, 1);
  }
  return restored;
}

std::vector<double> calculateResponses(
  const std::vector<TrackPure>& grid,
  const std::vector<Hit>& hits,
  double sharpness
)
{
  std::vector<double> responces(grid.size());
  for (unsigned int i = 0; i < grid.size(); ++i) 
  {
    auto& track = grid[i];
    for (const Hit& hit : hits)
    {
      responces[i] += exp(-getDistance(track, hit) / sharpness);
    }
  }
  return responces;
}

std::vector<Track> findHits(
  const std::vector<TrackPure>& tracks,
  const std::vector<Hit>& hits
)
{
  std::vector<Track> extendedTracks;
  extendedTracks.reserve(tracks.size());
  std::set<std::vector<uint32_t> > tracksSet;
  for (const TrackPure& track: tracks)
  {
    std::cerr << track.xOnZ0 << "," << track.yOnZ0 
      << "," << track.dxOverDz << "," 
      << track.dyOverDz << std::endl; 

    std::map<uint32_t, Hit> sensorsBest;
    for (const Hit& hit: hits) 
    {
      if (!sensorsBest.count(hit.sensorId) ||
          (getDistance(track, sensorsBest[hit.sensorId]) > 
          getDistance(track, hit)))
      {
        sensorsBest[hit.sensorId] = hit;
      }
    }
    Track extended;
    for (const auto& pair: sensorsBest)
    {
      if (pair.first < 0.6)
      {
        extended.addHit(pair.second.id);
      }
    }
    if (extended.hitsNum > 2)
    {
      extendedTracks.push_back(extended);
    }
  }
        
  return extendedTracks;
}

/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input  
 * @param output 
 */
int cpuRetinaInvocation(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output
)
{
  output.resize(input.size());
  for (size_t i = 0; i < input.size(); ++i)
  {
    EventInfo event = parseEvent(const_cast<const uint8_t*>(&(*input[i])[0]), input[i]->size());
    const std::vector<std::vector<double> > dimensions = generateDimensions(event);
    const std::vector<TrackPure> grid = generateGrid<TrackPure>(dimensions, generateTrackFromIndex);
    const std::vector<Hit>& hits = event.hits;
    const std::vector<double> responses = calculateResponses(grid, hits, RETINA_SHARPNESS_COEFFICIENT);
    const std::vector<TrackPure> restored = restoreTracks(dimensions, grid, responses);
    const std::vector<Track> tracksWithHits = findHits(restored, hits);
    auto answer = putTracksInOutputFormat(hits, tracksWithHits);
    output[i] = answer;
    //printSolution(tracksWithHits, hits, DEBUG);
  }
  return 0;
}
