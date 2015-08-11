#pragma once
#include<vector>
#include<functional>


std::vector<double> generateUniformDimension(double min, double max, size_t size);

int calculateGridSize(
  const std::vector<std::vector<double> >& dimensions, 
  const int noBoarders
);

void generateNextMultiIndex(
  std::vector<int>& multiIndex, 
  const std::vector<std::vector<double> >& dimensions, 
  int noBoarders
);


int multiIndexToIndex(
  const std::vector<int>& multiIndex, 
  const std::vector<std::vector<double> >& dimensions
);

std::vector<int> generateNeighboursIndexes(
  const std::vector<int>& multiIndex, 
  const std::vector<std::vector<double> >& dim
);

template<class T>
std::vector<T> generateGrid(
    const std::vector<std::vector<double> >& dimensions, 
    std::function<T(const std::vector<int>&, const std::vector<std::vector<double> >&)> generateOne
)
{
  std::vector<int> multiIndex(dimensions.size(), 0);
  int gridSize = calculateGridSize(dimensions, 0);
  std::vector<T> ans(gridSize);
  for (int i = 0; i < gridSize; ++i)
  {
    ans[i] = generateOne(multiIndex, dimensions);
    generateNextMultiIndex(multiIndex, dimensions, 0);
  }
  return ans;
}
