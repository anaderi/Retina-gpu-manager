// Automatically generated file
#pragma once
#include "FastTrackIpc/api/api.h"
#include <map>
#include <set>
#include <stdexcept>
#include <stdint.h>
#include <vector>
class IProtocol;
class Tracker {
public:
  // interface
  Tracker(IProtocol & protocol);
public:
  // service function wrappers
  void searchByPair(std::vector<Track> data, std::vector<int8_t> & result);
private:
  // data
  IProtocol & protocol;
};
