#include <stdexcept>
#include <algorithm>
#include <iostream>

#include "Grid.h"

std::vector<double> generateUniformDimension(double min, double max, size_t size)
{
  std::vector<double> dimension(size + 2);
  double quant = (max - min) / (size - 1);
  dimension[0] =  min - quant;
  for(size_t i = 0; i <= size; i++)
  {
    dimension[i + 1] = dimension[i] + quant;
  }
  return std::move(dimension);
}

std::vector<double> generateGaussDimension(double mean, double sd, size_t size)
{
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(mean, sd);
  std::vector<double> dimension(size);
  for (size_t i = 0; i < size; i++)
  {
    dimension[i] = distribution(generator);
  }
  std::sort(dimension.begin(), dimension.end());
  return std::move(dimension);
}

int calculateGridSize(
  const std::vector<std::vector<double> >& dimensions, 
  const int noBoarders
)
{
  int ans = 1;
  for (const std::vector<double>& dimension : dimensions)
    ans *= dimension.size() - 2 * noBoarders;
  return ans;
}

void generateNextMultiIndex(
  std::vector<int>& multiIndex, 
  const std::vector<std::vector<double> >& dimensions, 
  int noBoarders
)
{
  multiIndex[0]++;
  int i = 0;
  while (multiIndex[i] == dimensions[i].size() - noBoarders)
  {
    multiIndex[i] = noBoarders;
    ++multiIndex[i + 1];
    ++i;
  }
}

int multiIndexToIndex(
  const std::vector<int>& multiIndex, 
  const std::vector<std::vector<double> >& dimensions
)
{
  int index = multiIndex[multiIndex.size() - 1];
  for (int i = multiIndex.size() - 2; i >= 0; --i)
  {
    index = index * dimensions[i].size() + multiIndex[i];
  }
  if (index > calculateGridSize(dimensions, 0))
  {
    for (int i = 0; i < multiIndex.size(); ++i)
    {
      std::cerr << multiIndex[i] << " ";
    }
    std::cerr << std::endl;
  }
  return index;
}

std::vector<int> generateNeighboursIndexes(
  const std::vector<int>& multiIndex, 
  const std::vector<std::vector<double> >& dim
)
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
/*
int main()
{
  
}
*/
