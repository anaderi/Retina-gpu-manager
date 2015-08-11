#include <stdexcept>

#include "Grid.h"


double Dimension::getGridBoarder(int index) const
{
  if (index < 0 || index >= gridSize)
    throw std::logic_error("Bad index of gridBoarder");
  return min + quant * index;
}

int calculateGridSize(const std::vector<Dimension>& dimensions, const int noBoarders)
{
  int ans = 1;
  for (const Dimension& dimension : dimensions)
    ans *= dimension.gridSize - 2 * noBoarders;
  return ans;
}

void generateNextMultiIndex(std::vector<int>& multiIndex, const std::vector<Dimension>& dimensions, int noBoarders)
{
  multiIndex[0]++;
  int i = 0;
  while (multiIndex[i] == dimensions[i].gridSize - noBoarders)
  {
    multiIndex[i] = noBoarders;
    ++multiIndex[i + 1];
    ++i;
  }
}

int multiIndexToIndex(const std::vector<int>& multiIndex, const std::vector<Dimension>& dimensions)
{
  int index = multiIndex[multiIndex.size() - 1];
  for (int i = multiIndex.size() - 2; i >= 0; --i)
  {
    index = index * dimensions[i + 1].gridSize + multiIndex[i];
  }
  return index;
}

std::vector<int> generateNeighboursIndexes(const std::vector<int>& multiIndex, const std::vector<Dimension>& dim)
{
  std::vector<int> neighbours;
  neighbours.reserve(multiIndex.size() * 2);

  auto multiIndexCopy = multiIndex;
  for (int i = 0; i < multiIndex.size(); ++i)
  {
    for (int dx : {-1, 1})
    {
      multiIndexCopy[i] = multiIndex[i] + dx;
      neighbours.push_back(multiIndexToIndex(multiIndexCopy, dim));            
    }
    multiIndexCopy[i] = multiIndex[i];
  }
  return neighbours;
}

