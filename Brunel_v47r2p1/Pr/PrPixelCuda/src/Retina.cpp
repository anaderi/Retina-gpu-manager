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
  
  //const double threshold = 1.9;//getQuatile(responces, TRACK_THRESHOLD);
  for (int i = 0; i < gridSizeNoBoarders; ++i)
  {
    std::vector<int> neighbours = generateNeighboursIndexes(indexes, dimensions);
    int currentIndex = multiIndexToIndex(indexes, dimensions);
    bool isLocalMaximum = true;
    for (int neighbour : neighbours)
    {
      std::cerr << neighbour << " " << grid.size() << std::endl;
      if (responces[neighbour] > responces[currentIndex] + 1e-8)
      {
        isLocalMaximum = false;
      }
    }
    //isLocalMaximum = /*isLocalMaxumum ||*/ (responces[currentIndex] > threshold);
    if (isLocalMaximum && responces[currentIndex] > 1e-5)
    {
      TrackPure answer = grid[currentIndex] * responces[currentIndex];
      double sum_responce = responces[currentIndex];
      for (int neighbour : neighbours)
      {
        answer = answer + grid[neighbour] * responces[neighbour];
        sum_responce += responces[neighbour];
      }
      //std::cerr << responces[currentIndex] << "->" << sum_responce << std::endl;
      const TrackPure track = answer * (1.0 / sum_responce);
      restored.push_back(track);      
    /*
     *
     * 
     *  std::cerr << track.xOnZ0 << "," << track.yOnZ0 
        << "," << track.dxOverDz << "," 
        << track.dyOverDz << std::endl; 
     */
    }
    generateNextMultiIndex(indexes, dimensions, 1);
  }
  return restored;
}

std::vector<double> calculateResponses(
  const std::vector<TrackPure>& grid,
  const std::vector<Hit>& hits,
  const std::function<double(TrackPure,Hit)> distance,
  double sharpness
)
{
  std::vector<double> responces(grid.size());
  for (unsigned int i = 0; i < grid.size(); ++i) 
  {
    auto& track = grid[i];
    for (const Hit& hit : hits)
    {
      responces[i] += exp(-distance(track, hit) / sharpness);
    }
  }
  return responces;
}
std::vector<std::pair<double, Hit> > findHitsOnSensor(  
  const std::vector<Hit>& hits,
  const TrackPure& track,
  const std::function<double(TrackPure,Hit)> distance
)
{
  std::vector<std::pair<double, Hit> > distances;
  std::map<uint32_t, Hit> sensorsBest;

  for (size_t i = 0; i < hits.size(); ++i)
  {
    {
      const Hit& hit = hits[i];
      if (!sensorsBest.count(hit.sensorId) ||
          (distance(track, sensorsBest[hit.sensorId]) > 
          distance(track, hit)))
      {
        sensorsBest[hit.sensorId] = hit;
      }
    }
  }
  for (const auto& pair: sensorsBest)
  {
    distances.emplace_back(getDistance(track, pair.second), pair.second);
  }
  return distances;
}

std::vector<Track> findHits(
  const std::vector<TrackPure>& tracks,
  const std::vector<Hit>& hits
)
{
  std::vector<Track> extendedTracks;
  extendedTracks.reserve(tracks.size());
  std::set<std::vector<uint32_t> > tracksSet;
  std::vector<bool> used(hits.size(), false);
  int cnt = 0;
  for (const TrackPure& track: tracks)
  {
    if (++cnt % 1000 == 0)
    {
      std::cerr << cnt << " tracks proceed" << std::endl;
    }
    std::vector<std::pair<double, Hit> > distances = findHitsOnSensor(hits, track, getDistance) ;
    Track extended;
    std::sort(distances.begin(), distances.end(), 
      [](const std::pair<double, Hit>& a, const std::pair<double, Hit>& b) -> bool
      {
        return a.first < b.first;
      }
    );
    {
      double zStart = distances[0].second.z;
      extended.addHit(distances[0].second.id);
      for (size_t i = 1; i < 3; /*distances.size()*/ i++)
      {
        if (isFit(track, distances[i].second, zStart))
        {
          extended.addHit(distances[i].second.id);
        }
      }
    }
    if (extended.hitsNum > 2)
    {
      /*for (size_t i = 0; i < extended.hitsNum; ++i)
      {
        used[extended.hits[i]] = true;
      }*/
      extendedTracks.push_back(extended);
    }
  }
        
  return extendedTracks;
}
/*
TrackPure makeLocalMaximum(TrackPure track, const std::vector<Hit>& hits, double rs) {
  std::vector<TrackPure> steps = {
    TrackPure(1, 0, 0, 0),
    TrackPure(0, 1, 0, 0),
    TrackPure(0, 0, 1, 0),
    TrackPure(0, 0, 0, 1)
  };
  double len = 1e-3;
  while (len > 1e-7)
  {
    bool update = false;
    TrackPure best = track;
    double bestResponce = calculateResponses({track}, hits, rs)[0];
    for (double dx: {-len, len})
    {
      for (TrackPure step: steps)
      {
        double current = calculateResponses({track + step * dx}, hits, rs)[0];
        if ( current > bestResponce + 1e-7)
        {
          update = true;
          best = track + step * dx;
          bestResponce = current;
        }
      }
    }
    if (!update)
    {
      len /= 2;
    }
    else
    {
      track = best;
      len *= 2;
    }
    //std::cerr << len << std::endl;
    //std::cerr << bestResponce << std::endl;
  }
  return track;
}
 * */

void outputBest(
  const std::vector<TrackPure>& tracks,
  const std::vector<Hit>& hits,
  const std::function<double(TrackPure,Hit)> distance,
  std::string name
)
{
  std::ofstream myfile;
  myfile.open (name);
  for (const TrackPure& track: tracks)
  {
    auto bst = findHitsOnSensor(hits, track, distance);
    for (auto p: bst)
    {
      myfile << p.second.id << " ";
    }
    myfile << std::endl;
  }
  myfile.close();
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
    const std::vector<std::vector<double> > dimDx = {
      generateUniformDimension(-1, 1, 198),
      generateUniformDimension(0, 0, 1),
      generateUniformDimension(-0.3, 0.3, 198),
      generateUniformDimension(0, 0, 1)
    };
    const std::vector<std::vector<double> > dimDy = {
      generateUniformDimension(0, 0, 1),
      generateUniformDimension(-1, 1, 198),
      generateUniformDimension(0, 0, 1),
      generateUniformDimension(-0.3, 0.3, 198)
    };

    std::cerr << "Gx" << std::endl;
    const std::vector<TrackPure> gridDx = generateGrid<TrackPure>(dimDx, generateTrackFromIndex);
    std::cerr << "Gy" << std::endl;
    const std::vector<TrackPure> gridDy = generateGrid<TrackPure>(dimDy , generateTrackFromIndex);

    const std::vector<Hit>& hits = event.hits;
    std::cerr << "calc" << std::endl;
    const std::vector<double> responsesDx = calculateResponses(gridDx, hits, getDistanceDx, RETINA_SHARPNESS_COEFFICIENT);
    std::cerr << "restore" << std::endl;
    const std::vector<TrackPure> restoredDx = restoreTracks(dimDx, gridDx, responsesDx);

    std::cerr << "calc" << std::endl;
    const std::vector<double> responsesDy = calculateResponses(gridDy, hits, getDistanceDy, RETINA_SHARPNESS_COEFFICIENT);
    std::cerr << "restore" << std::endl;
    const std::vector<TrackPure> restoredDy = restoreTracks(dimDy, gridDy, responsesDy);
    std::cerr << "output" << std::endl;
    {
      outputBest(restoredDx, hits, getDistanceDx, "dx.txt");
      outputBest(restoredDy, hits, getDistanceDy, "dy.txt");      
    }
    std::cerr << "in Find" << std::endl;
    std::vector<TrackPure> tracks;
    for (const TrackPure& dx : restoredDx)
    {
      std::vector<std::pair<double, Hit> > xHits = findHitsOnSensor(hits, dx, getDistance) ;
      for (const TrackPure& dy: restoredDy)
      {
        std::vector<std::pair<double, Hit> > yHits = findHitsOnSensor(hits, dx, getDistance) ;
        int intersection = 0;

        for (int j = 0; j < xHits.size(); j++)
        {
          if (xHits[j].first + yHits[j].first < PARAM_TOLERANCE && xHits[j].second.id == yHits[j].second.id)
          {
            intersection++;
          }
        }
        if (intersection > 2)
        {
          tracks.push_back(dx + dy);
        }
      }
    }
    auto answer = putTracksInOutputFormat(hits, findHits(tracks, hits));
    output[i] = answer;
    //printSolution(tracksWithHits, hits, DEBUG);
  }
  return 0;
}
