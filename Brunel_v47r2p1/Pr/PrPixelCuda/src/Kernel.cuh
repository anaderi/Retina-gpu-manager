#pragma once

#ifdef CUDA_TEST
#include "cuda-cleanup.h"
#endif

#include "Definitions.h"

__device__ float fitHits(Hit& h0, Hit& h1, Sensor& s0, Sensor& s1);
__device__ float fitHitToTrack(Track& t, Hit& h1, Sensor& s1);
__device__ void acceptTrack(Track& t, TrackFit& fit, Hit& h0, Hit& h1, Sensor& s0, Sensor& s1, int h0_num, int h1_num);
__device__ void updateTrack(Track& t, TrackFit& fit, Hit& h1, Sensor& s1, int h1_num);
__device__ void updateTrackCoords(Track& t, TrackFit& fit);
__device__ float trackChi2(Track& t);
__device__ float hitChi2(Track& t, Hit& h, int hit_z);

__global__ void prepareData(char* input, int* _prevs, int* _nexts);
__global__ void gpuKalman(Track* tracks, bool* track_holders);
__global__ void postProcess(Track* tracks, bool* track_holders, int* track_indexes, int* num_tracks, int* tracks_to_process);

__global__ void gpuKalmanBalanced(Span  * spans, Fit * fittings);
__global__ void consolidateHits(Fit * fittings, int n, Track * tracks, bool * trackHolders);

__device__ void addTrack(const Fit & fit, Track * tracks, bool * trackHolders);
