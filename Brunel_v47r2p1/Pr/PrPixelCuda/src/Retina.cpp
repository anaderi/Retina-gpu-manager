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
  const std::vector<Dimension>& dim
)
{
  return TrackPure(
    dim[0].getGridBoarder(multiIndex[0]),
    dim[1].getGridBoarder(multiIndex[1]),
    dim[2].getGridBoarder(multiIndex[2]),
    dim[3].getGridBoarder(multiIndex[3])
  );
}

std::vector<TrackPure> restoreTracks(
  const std::vector<Dimension>& dimensions,
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
/*std::vector<TrackPure> climbingTheHill(std::vector<TrackPure> tracks)
{
  const std::vector<Dimension> dimensions =
    {
      Dimension(-1, 1, 3),
      Dimension(-1, 1, 3),
      Dimension(-1, 1, 3),
      Dimension(-1, 1, 3)
    };
  const std::vector<TrackPure> grid = generateGrid<TrackPure>(dimensions, generateTrackFromIndex);

  double step = 1;
  while (step > 1e-4) 
  {
    for (int i = 0; i < tracks.size(); i++)
    {
      TrackPure best = tracks[i];
      double bestResponce = calculateResponses
    }
  }
}
*/
std::vector<Track> findHits(
  const std::vector<TrackPure>& tracks,
  const std::vector<Hit>& hits
)
{
  std::vector<Track> extendedTracks;
  extendedTracks.reserve(tracks.size());
  std::set<std::vector<uint32_t> > tracksSet;
  for (const TrackPure& track1: tracks)
  {
    TrackPure track = track1;
    for (int j = 1; j < 2; j++)
    {
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
      double sumDistance = 0;
      std::vector<std::pair<double, Hit> > distances;
      for (const auto& pair: sensorsBest)
      {
        distances.emplace_back(getDistance(track, pair.second), pair.second);
      }
      sort(distances.begin(), distances.end(), 
        [](const std::pair<double, Hit>& a, const std::pair<double, Hit>& b) -> bool
        {
          return a.first < b.first;
        });
      if (j == 0)
      {
        Hit a = distances[0].second;
        Hit b = distances[1].second;
        float tx = (b.x - a.x) / (b.z - a.z);
        float ty = (b.y - a.y) / (b.z - a.z);
        float x0 = a.x - a.z * tx;
        float y0 = a.y - a.z * ty;
        
        track = TrackPure(x0, y0, tx, ty);
      }
      else
      {
        Track extended;
        for (int i = 0; i < 3; i++)
        {
          extended.addHit(distances[i].second.id);
          //extended.addHit(rand() % hits.size());
          //std::cerr << "putting random" << std::endl;
        }
        extendedTracks.push_back(extended);
      }
      
//      {
//        update
//      }
//      if (!tracksSet.count(trackVc))
//      {
//        //std::cerr << "hits" << extended.hitsNum << "average distance from hits: " 
//        //  << sumDistance / extended.hitsNum << std::endl;
//        tracksSet.insert(trackVc);
//        extendedTracks.push_back(extended);
//      }
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
    const std::vector<Dimension> dimensions = generateDimensions(event);
    const std::vector<TrackPure> grid = generateGrid<TrackPure>(dimensions, generateTrackFromIndex);
    const std::vector<Hit>& hits = event.hits;
    const std::vector<double> responses = calculateResponses(grid, hits, RETINA_SHARPNESS_COEFFICIENT);
    const std::vector<TrackPure> restored = restoreTracks(dimensions, grid, responses);
    const std::vector<Track> tracksWithHits = findHits(restored, hits);
    std::cerr << "T Size = " << tracksWithHits.size() << std::endl;
    auto answer = putTracksInOutputFormat(hits, tracksWithHits);
    output[i] = answer;
    std::cerr << "Size = " << answer.size() << std::endl;
    //printSolution(tracksWithHits, hits, DEBUG);
  }
  return 0;
}
