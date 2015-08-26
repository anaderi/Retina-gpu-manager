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
  
  const double threshold = 1.9;//getQuatile(responces, TRACK_THRESHOLD);
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
  std::vector<bool> used(hits.size(), false);
  for (const TrackPure& track: tracks)
  {

    std::map<uint32_t, Hit> sensorsBest;
    for (size_t i = 0; i < hits.size(); ++i)
    {
      //if (!used[i])
      {
        const Hit& hit = hits[i];
        if (!sensorsBest.count(hit.sensorId) ||
            (getDistance(track, sensorsBest[hit.sensorId]) > 
            getDistance(track, hit)))
        {
          sensorsBest[hit.sensorId] = hit;
        }
      }
    }
    std::vector<std::pair<double, Hit> > distances;
    Track extended;
    for (const auto& pair: sensorsBest)
    {
      distances.emplace_back(getDistance(track, pair.second), pair.second);
    }
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
        //if (isFit(track, distances[i].second, zStart))
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
    //const std::vector<double> responses = calculateResponses(grid, hits, RETINA_SHARPNESS_COEFFICIENT);
    //const std::vector<TrackPure> restored = restoreTracks(dimensions, grid, responses);
    const std::vector<TrackPure> restored;
    
    /*{
      std::ofstream myfile;
      myfile.open ("tracks.csv");
      myfile << "x0,y0,dx,dy,reps" << std::endl;
      for (const TrackPure& track: restored)
      {
        myfile << track.xOnZ0 << "," << track.yOnZ0 << "," 
          << track.dxOverDz << "," << track.dyOverDz << ","
          << calculateResponses({track}, hits, RETINA_SHARPNESS_COEFFICIENT)[0]<< std::endl;
      }
      myfile.close();
    }*/
    std::vector<TrackPure> maxims;

    {
      std::ofstream myfile;
      myfile.open ("maxims.csv");
      myfile << "x0,y0,dx,dy,reps" << std::endl;
      std::default_random_engine generator;
      std::normal_distribution<double> dx(0, 0.3);
      std::normal_distribution<double> x0(0, 10);
      std::uniform_int_distribution<int> distribution(0, hits.size());  

      int cnt = 0, nw = 0;
      while (cnt < 10000 && maxims.size() < 10000)
      //for (size_t i = 0 ; i < )
      {
        TrackPure track;
        {
          track = makeLocalMaximum(
            TrackPure(
              x0(generator),
              x0(generator),
              dx(generator),
              dx(generator)
            ),
            hits, 
            RETINA_SHARPNESS_COEFFICIENT
            );
        }
        if (std::all_of(maxims.begin(), maxims.end(), 
          [&](TrackPure& old)
          {
            if (fabs(old.xOnZ0 - track.xOnZ0) > 0.1)
              return true;
            if (fabs(old.yOnZ0 - track.yOnZ0) > 0.1)
              return true;
            if (fabs(old.dxOverDz - track.dxOverDz) > 0.01)
              return true;
            if (fabs(old.dyOverDz - track.dyOverDz) > 0.01)
              return true;
            return false;
          }) 
          )
        {
          maxims.push_back(track);
          nw = 0;
          myfile << track.xOnZ0 << "," << track.yOnZ0 << "," 
            << track.dxOverDz << "," << track.dyOverDz << ","
            << calculateResponses({track}, hits, RETINA_SHARPNESS_COEFFICIENT)[0]<< std::endl;
        }
        else
        {
          nw++;
        }
        if (++cnt % 10 == 0)
        {
          std::cerr << cnt << " tracks procceed:" << maxims.size() <<  std::endl;
        }
       }
      
    
      myfile.close();
    }
    const std::vector<double> resp = calculateResponses(maxims, hits, RETINA_SHARPNESS_COEFFICIENT);
    std::vector<TrackPure> final;
    {
      double threshold = -1;// getQuatile(resp, 0);
      for(size_t i = 0; i < maxims.size(); i++)
        if (resp[i] > threshold)
          final.push_back(maxims[i]);
    }
    const std::vector<Track> tracksWithHits = findHits(final, hits);
    auto answer = putTracksInOutputFormat(hits, tracksWithHits);
    output[i] = answer;
    //printSolution(tracksWithHits, hits, DEBUG);
  }
  return 0;
}
