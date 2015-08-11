#pragma once
#include<vector>
#include<functional>

struct Dimension
{
  double min;
  double max;
  double quant;
  size_t gridSize;
  Dimension(double min_, double max_, size_t gridSize_)
  {
    if (gridSize_ < 2)
      throw std::runtime_error("Grid size must be > 1, given:" + gridSize_);
    quant = (max_ - min_) / (gridSize_ - 1);
    min = min_ - quant;
    max = max_ + quant;
    gridSize = gridSize_ + 2;
  }
  
  double getGridBoarder(int index) const;
};

int calculateGridSize(const std::vector<Dimension>& dimensions, const int noBoarders);

void generateNextMultiIndex(std::vector<int>& indexes, const std::vector<Dimension>& dim, int noBoarders);

int multiIndexToIndex(const std::vector<int>& indexes, const std::vector<Dimension>& dim);

std::vector<int> generateNeighboursIndexes(const std::vector<int>& indexes, const std::vector<Dimension>& dim);

template<class T>
std::vector<T> generateGrid(
    const std::vector<Dimension>& dimensions, 
    std::function<T(const std::vector<int>&, const std::vector<Dimension>&)> generateOne
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
