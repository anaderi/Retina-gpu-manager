#include "Physics.h"

#include<iostream>
#include <limits>
TrackPure operator*(const TrackPure& one, const double alpha)
{
  return TrackPure(
      one.xOnZ0 * alpha,
      one.yOnZ0 * alpha,
      one.dxOverDz * alpha,
      one.dyOverDz * alpha
  );
}

TrackPure operator+(const TrackPure& one, const TrackPure& other) 
{
  return TrackPure(
    one.xOnZ0 + other.xOnZ0,
    one.yOnZ0 + other.yOnZ0,
    one.dxOverDz + other.dxOverDz,
    one.dyOverDz + other.dyOverDz
  );
}
TrackPure::TrackPure(const Hit& a, const Hit& b)
{
  dxOverDz = (b.x - a.x) / (b.z - a.z);
  dyOverDz = (b.y - a.y) / (b.z - a.z);
  xOnZ0 = a.x - a.z * dxOverDz;
  yOnZ0 = a.y - a.z * dyOverDz;
}

inline double square(double x)
{
  return x * x;
}

double getDistance(const TrackPure& track, const Hit& hit) noexcept
{
  return square(hit.x - track.xOnZ0 - track.dxOverDz * hit.z) +
         square(hit.y - track.yOnZ0 - track.dyOverDz * hit.z);
}

std::vector<std::vector<double> > generateDimensions(const EventInfo& event)
{
  double minX0, maxX0;
  double minY0, maxY0;
  double minDx, maxDx;
  double minDy, maxDy;
  minX0 = minY0 = minDx = minDy = std::numeric_limits<double>::max();
  maxX0 = maxY0 = maxDx = maxDy = std::numeric_limits<double>::min();
  
  for (size_t i = 0; i < event.hits.size(); ++i)
  {
    for (size_t j = 0; j < i; j++)
    {
      if (event.hits[i].sensorId != event.hits[j].sensorId)
      {
        TrackPure t(event.hits[i], event.hits[j]);
        minX0 = std::min(minX0, t.xOnZ0);
        maxX0 = std::max(maxX0, t.xOnZ0);
        minY0 = std::min(minY0, t.yOnZ0);
        maxY0 = std::max(maxY0, t.yOnZ0);
        minDx = std::min(minDx, t.dxOverDz);
        maxDx = std::max(maxDx, t.dxOverDz);
        minDy = std::min(minDy, t.dyOverDz);
        maxDy = std::max(maxDy, t.dyOverDz);
      }
    }
  }
  std::cerr << "Dimensions generated" << std::endl;
  return std::vector<std::vector<double> > {
        generateUniformDimension(minX0, maxX0, GRID_SIZE_X_ON_Z0),
        generateUniformDimension(minY0, maxY0, GRID_SIZE_Y_ON_Z0),
        generateUniformDimension(minDx, maxDx, GRID_SIZE_DX_OVER_DZ),
        generateUniformDimension(minDy, maxDy, GRID_SIZE_DY_OVER_DZ)
  };
}
