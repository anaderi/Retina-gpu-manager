#include "Physics.h"

#include <iostream>
#include <limits>
#include <cmath>
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
  /*return std::vector<std::vector<double> > {
        generateUniformDimension(-1, 1, GRID_SIZE_X_ON_Z0),
        generateUniformDimension(-1, 1, GRID_SIZE_Y_ON_Z0),
        generateUniformDimension(-0.01, 0.01, GRID_SIZE_DX_OVER_DZ),
        generateUniformDimension(-0.01, 0.01, GRID_SIZE_DY_OVER_DZ)
  };*/
  return std::vector<std::vector<double> > {
        generateGaussDimension(0, 12, GRID_SIZE_X_ON_Z0),
        generateGaussDimension(0, 12, GRID_SIZE_Y_ON_Z0),
        generateGaussDimension(0, 0.1, GRID_SIZE_DX_OVER_DZ),
        generateGaussDimension(0, 0.1, GRID_SIZE_DY_OVER_DZ)
  };
  
}

bool isFit(const TrackPure& track, const Hit& hit, double zStart) noexcept
{
  double dx = fabs(hit.x - track.xOnZ0 - track.dxOverDz * hit.z);
  double dy = fabs(hit.y - track.yOnZ0 - track.dyOverDz * hit.z);
  double scatterNum = dx * dx + dy * dy;
  double scatterNorm = 1 / (hit.z - zStart);
  double scatter = scatterNum * scatterNorm * scatterNorm;
  return dx < PARAM_TOLERANCE && dy < PARAM_TOLERANCE && scatter < SCATTER_TOLERANCE;
}