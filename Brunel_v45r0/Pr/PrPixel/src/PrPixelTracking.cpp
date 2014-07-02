// Gaudi
#include "GaudiKernel/AlgFactory.h"
// LHCb
#include "Event/Track.h"
#include "Event/StateParameters.h"
// Local
#include "PrPixelTracking.h"
#include "PrPixelSerialization.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrPixelTracking
//
// 2011-12-16 : Olivier Callot
//-----------------------------------------------------------------------------

DECLARE_ALGORITHM_FACTORY(PrPixelTracking)

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrPixelTracking::PrPixelTracking(const std::string& name,
                                 ISvcLocator* pSvcLocator) :
#ifdef DEBUG_HISTO
    GaudiTupleAlg(name, pSvcLocator),
#else
    GaudiAlgorithm(name, pSvcLocator),
#endif
    m_debugTool(NULL) { 

  declareProperty("OutputTracksName", m_outputLocation = LHCb::TrackLocation::Velo);
  declareProperty("MaxXSlope", m_maxXSlope = 0.400);
  declareProperty("MaxYSlope", m_maxYSlope = 0.400);
  // Tolerance window when adding hits
  declareProperty("ExtraTol", m_extraTol = 0.6 * Gaudi::Units::mm);
  // Number of modules without a hit after which to stop extrapolation
  declareProperty("MaxMissed", m_maxMissed = 3); 

  // Acceptance criteria for adding new hits
  declareProperty("MaxScatter", m_maxScatter = 0.004);
  
  // Acceptance criteria for track candidates
  // Max. chi2 for 3-hit tracks
  declareProperty("MaxChi2Short", m_maxChi2Short = 20.0); 
  // Min. fraction of unused hits
  declareProperty("FractionUnused", m_fractionUnused = 0.5);

  // Flag to clear hits (for rerunning in same event) 
  declareProperty("ClearHits", m_clearHits = false);  

  declareProperty("UseSlopeCorrection", m_useSlopeCorrection = false);

  // Parameters for debugging
  declareProperty("DebugToolName", m_debugToolName = "");
  declareProperty("WantedKey", m_wantedKey = -100);
  declareProperty("TimingMeasurement", m_doTiming = false);

  // Parameters for Kalman fit
  declareProperty("ClosestToBeamStateKalmanFit", m_stateClosestToBeamKalmanFit = true);
  declareProperty("EndVeloStateKalmanFit", m_stateEndVeloKalmanFit = false);
  declareProperty("AddFirstLastMeasurementStatesKalmanFit", m_addStateFirstLastMeasurementKalmanFit = false);

}

//=============================================================================
// Destructor
//=============================================================================
PrPixelTracking::~PrPixelTracking() {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrPixelTracking::initialize() {

  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;
  m_isDebug = msgLevel(MSG::DEBUG);
  // Setup the hit manager.
  m_hitManager = tool<PrPixelHitManager>("PrPixelHitManager");
  m_hitManager->useSlopeCorrection(m_useSlopeCorrection);
  // Setup the debug tool.
  if ("" != m_debugToolName) m_debugTool = tool<IPatDebugTool>(m_debugToolName);
  // Setup the timing measurement.
  if (m_doTiming) {
    m_timerTool = tool<ISequencerTimerTool>("SequencerTimerTool/Timer", this);
    m_timeTotal = m_timerTool->addTimer("Total");
    m_timerTool->increaseIndent();
    m_timePrepare = m_timerTool->addTimer("Prepare");
    m_timePairs = m_timerTool->addTimer("Find by pairs");
    m_timeFinal = m_timerTool->addTimer("Store tracks");
    m_timerTool->decreaseIndent();
  }
#ifdef DEBUG_HISTO
  setHistoTopDir("VP/");
#endif

  gpuService = svc<IGpuService>("GpuService", true);

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================

/// Callback function used with GpuService::submitData.
/// allocTracks takes the size of received data and a pointer to a GpuTrack
/// vector. The received data is assumed to consist of an array of GpuTrack
/// objects. allocTracks reserves enough space to store the received tracks and
/// returns a pointer to the start of that memory.
void * allocTracks(size_t size, void * param)
{
  typedef std::vector<uint8_t> Data;
  Data & tracks = *reinterpret_cast<Data*>(param);
  tracks.resize(size);
  return &tracks[0]; // size is strictly positive
}

StatusCode PrPixelTracking::execute() {
  // Code with minimal instructions to set it off to the GPU
  // and perform check afterwards.
  // Also timing unaffected.
  
  if (m_doTiming) {
    m_timerTool->start(m_timeTotal);
    m_timerTool->start(m_timePrepare);
  }

  // Clean and build hits from clusters
  // From what PrPixel has:
  //  PrPixel format -> SoA (done alongside buildHits) -> AoS
  if (m_clearHits) m_hitManager->clearHits();
  m_hitManager->buildHits();

  // Do some typecasting into a format understandable by GPU
  m_hitManager->m_serializer.serializeEvent(m_serializedEvent);
  
  if (m_serializedEvent.empty())
    info() << "--- Serialized event is empty! This should not happen!" << endmsg;
  
  if (m_isDebug){
    info() << "--- Submitting data to gpuService" << endmsg;
    info() << "--- serializedEvent: 0x" << std::hex << (long long)&m_serializedEvent[0] <<
      std::dec << ", size " << m_serializedEvent.size() << endmsg;
  }

  if (m_doTiming) m_timerTool->stop(m_timePrepare);

  if (m_doTiming) m_timerTool->start(m_timePairs);

  // Perform search on the GPU
  std::vector<uint8_t> trackCollection;

  try {
    gpuService->submitData("PrPixelCudaHandler", &m_serializedEvent[0], m_serializedEvent.size(), allocTracks, &trackCollection);

    // There is no need to deserialize the tracks,
    // and then convert them to "CPU" tracks. These two
    // can be done in one go (that's what "deserializeTracks" should do, and
    // it's PrPixel dependent)
    // if (m_isDebug)
    //   info() << "--- Deserializing gpu tracks" << endmsg;


    // Conversion from track to PrPixelTrack
    m_hitManager->m_serializer.deserializeTracks(trackCollection, m_tracks);

    // Additional minor things, which should be left out afterwards
    // for (PrPixelTracks::iterator it = m_tracks.begin(); it != m_tracks.end(); it++){
    //   if ( it->hits().size() > 3 )
    //     it->tagUsedHits();
    // }
  } catch (const std::exception & e) {
    error() << "submission failed; " << e.what() << std::endl;
  } catch (...) {
    error() << "submission failed; reason unknown" << std::endl;
  }

  if (m_doTiming) m_timerTool->stop(m_timePairs);

  // Convert temporary tracks to LHCb tracks.
  if (m_doTiming) m_timerTool->start(m_timeFinal);
  makeLHCbTracks();
  if (m_doTiming) {
    m_timerTool->stop(m_timeFinal);
    m_timerTool->stop(m_timeTotal);
  }


#ifdef DEBUG_HISTO
  for (unsigned int i = m_hitManager->firstModule(); i < m_hitManager->lastModule(); ++i) {
    PrPixelHits::iterator ith;
    for (ith = m_hitManager->hits(i).begin(); ith != m_hitManager->hits(i).end(); ++ith) {
      PrPixelHit* iC = *ith;
      const double x = iC->x();
      const double y = iC->y();
      const double z = iC->z();
      const double r = sqrt(x * x + y * y);
      if (!iC->isUsed()) {
        plot3D(x, y, z, "UnusedHits3D", "Distribution of UnusedHits",-50.0,50.0,-50.0,50.0,-500.0,800.0,100,100,200);
        plot2D(r, z, "UnusedHits_rz", "Distribution of Unused Hits", 0.0, 60.0,-500.0,800.0,100,100);
        plot2D(x, y, "UnusedHits2D", "Distribution of Unused Hits",-50.0,50.0,-50.0,50.0,100,100);
      }
      plot3D(x, y, z, "Hits_3D", "3D Distribution of Hits",-50.0,50.0,-50.0,50.0,-500.0,800.0,100,100,200);
      plot2D(r, z, "Hits_RZ", "RZ Distribution of Hits", 0.0, 60.0,-500.0,800.0,100,100);
      plot2D(x, y, "Hits_2D", "2D Distribution of Hits",-50.0,50.0,-50.0,50.0,100,100);
    }
  }

  const unsigned int nbHits = m_hitManager->nbHits();
  plot(nbHits, "HitsPerEvent", "Number of hits per event", 0.0, 8000.0, 80);
  plot(m_tracks.size(), "TracksPerEvent", "Number of tracks per event", 0.0,  800.0, 80);
  if (nbHits > 0) {
    const unsigned int nbHitsUsed = m_hitManager->nbHitsUsed();
    plot((100.0*nbHitsUsed)/nbHits, "PercentUsedHitsPerEvent", "Percent of hits assigned to tracks", 0.0, 100.0, 100);
  }
  for (PrPixelTracks::const_iterator itT = m_tracks.begin(); itT != m_tracks.end(); ++itT) {
    if ((*itT).size() <= 3) continue;
    // Calculate radius at first and last hit (assume that hits are sorted by module)
    const double x1 = (*itT).hits.front()->x();
    const double y1 = (*itT).hits.front()->y();
    const double r1 = sqrt(x1 * x1 + y1 * y1);
    const double x2 = (*itT).hits.back()->x();
    const double y2 = (*itT).hits.back()->y();
    const double r2 = sqrt(x2 * x2 + y2 * y2);
    const double minR = r1 > r2 ? r2 : r1;
    const double maxR = r1 > r2 ? r1 : r2;
    plot(minR, "MinHitRadiusPerTrack", "Smallest hit radius [mm] per track (of 4 or more hits)", 0.0, 50.0, 100);
    plot(maxR, "MaxHitRadiusPerTrack", "Largest hit radius [mm] per track (of 4 or more hits)",  0.0, 50.0, 100);
  }
#endif

  return StatusCode::SUCCESS;
}

//=============================================================================
// Extend track towards smaller z, 
// on both sides of the detector as soon as one hit is missed.
//=============================================================================
void PrPixelTracking::extendTrack(const PrPixelHit* h1, 
                                  const PrPixelHit* h2) {

  // Initially scan every second module (stay on the same side).
  int step = 2;
  // Start two modules behind the last one.
  int next = h2->module() - step; 
  // Count modules without hits found.
  unsigned int nbMissed = 0;
  while (next >= 0) {
    PrPixelModule* module = m_hitManager->module(next);
    // Attempt to add new hits from this module with given tolerances
    PrPixelHit* h3 = bestHit(module, m_extraTol, m_maxScatter, h1, h2);
    if (h3) {
      m_track.addHit(h3);
      // Reset missed hit counter.
      nbMissed = 0;
      // Update the pair of hits to be used for extrapolating.
      h1 = h2;
      h2 = h3;
    } else {
      // No hits found.
      if (step == 2) {
        // Look on the other side.
        module = m_hitManager->module(next + 1);
        h3 = bestHit(module, m_extraTol, m_maxScatter, h1, h2);
        if (!h3) {
          nbMissed += step;
        } else {
          m_track.addHit(h3);
          h1 = h2;
          h2 = h3;
        }
        // Switch to scanning every module (left and right). 
        step = 1; 
      } else {
        ++nbMissed;
      }
    }
    if (m_maxMissed < nbMissed) break;
    next -= step;
  }

}

//=========================================================================
//  Search starting with a pair of consecutive modules.
//=========================================================================
void PrPixelTracking::searchByPair() {

  // Get the range of modules to start search on,
  // starting with the one at largest Z. 
  const int lastModule  = m_hitManager->lastModule();
  const int firstModule = m_hitManager->firstModule() + 2;
  for (int sens0 = lastModule; firstModule <= sens0; sens0 -= 1) { 
    // Pick-up the "paired" module one station backwards
    const int sens1 = sens0 - 2;
    PrPixelModule* module0 = m_hitManager->module(sens0);
    PrPixelModule* module1 = m_hitManager->module(sens1);
    double z0 = module0->z();
    double z1 = module1->z();
    double dz = z0 - z1;
    // Does it make sense to skip pairs very far apart ?
    // if( fabs(dz) > 60.0) continue;                                      
#ifdef DEBUG_HISTO
    plot(dz, "SeedPairDeltaZ", "Separation in Z [mm] between the seed pair modules", -200.0, +200.0, 400);
#endif
    // Calculate the search window from the slope limits. 
    const double dxMax = m_maxXSlope * fabs(dz);
    const double dyMax = m_maxYSlope * fabs(dz);
    // Loop over hits in the first module (larger Z) in the pair.
    PrPixelHits::const_iterator ith0;
    PrPixelHits::const_iterator last0 = module0->hits().end();
    PrPixelHits::const_iterator first1 = module1->hits().begin();
    for (ith0 = module0->hits().begin(); last0 != ith0; ++ith0) {
      // Skip hits already assigned to tracks.
      if ((*ith0)->isUsed()) continue;                                            
      const double x0 = (*ith0)->x();
      const double y0 = (*ith0)->y();
      // Calculate x-pos. limits on the other module.
      const double xMin = x0 - dxMax;
      const double xMax = x0 + dxMax;
      if (m_debugTool && matchKey(*ith0)) {
        info() << format("s1%3d xMin%9.3f xMax%9.3f ", sens1, xMin, xMax);
        printHit(*ith0, "St0");
      }
      // Loop over hits in the second module (smaller Z) in the pair.
      PrPixelHits::const_iterator ith1;
      PrPixelHits::const_iterator last1 = module1->hits().end();
      for (ith1 = first1; last1 != ith1; ++ith1) { 
        const double x1 = (*ith1)->x();
        // Skip hits below the X-pos. limit.
        if (x1 < xMin) {
          first1 = ith1 + 1;
          continue;
        } 
        // Stop search when above the X-pos. limit.
        if (x1 > xMax) break;
        // Skip hits already assigned to tracks.
        if ((*ith1)->isUsed()) continue;
        // Check y compatibility.
        const double y1 = (*ith1)->y();
        // Skip hits out of Y-pos. limit.
        if (fabs(y1 - y0) > dyMax) continue;

        m_debug = m_isDebug;
        if (m_debugTool) {
          if (matchKey(*ith0) && matchKey(*ith1)) m_debug = true;
          if (m_debug) {
            info() << format("s1%3d dxRel %7.3f dyRel %7.3f    ", 
                             sens1, (x1-xMin)/(xMax-xMin), fabs((*ith1)->y()-y0)/dyMax);
            printHit(*ith1);
          }
        }
        // Make a seed track out of these two hits.
        m_track.set(*ith0, *ith1);
        // Extend the seed track towards smaller Z.
        extendTrack(*ith0, *ith1);
        if (m_track.hits().size() < 3) continue;

        // Final checks
        if (m_track.hits().size() == 3) {
          // If only 3 hits, all should be unused and chi2 good.
          if (m_track.nbUnused() != 3) {
            if (m_debug) {
              info() << "  -- reject, only " << m_track.nbUnused() << " unused hits." << endmsg;
              printTrack(m_track); 
            }
            continue;
          }
          if (m_track.chi2() > m_maxChi2Short) {
            if (m_debug) {
              info() << " -- reject, chi2 " << m_track.chi2() << " too high." << endmsg;
              printTrack(m_track);
            }
            continue;
          }
        } else {
          if (m_track.nbUnused() < m_fractionUnused * m_track.hits().size()) {
            if (m_debug) {
              info() << "  -- reject, only " << m_track.nbUnused() << "/" 
                     << m_track.hits().size() << " hits are unused." << endmsg;
              printTrack(m_track); 
            }
            continue;
          }
        }

        m_tracks.push_back(m_track);
        if (m_debug) {
          info() << "=== Store track Nb " << m_tracks.size() << endmsg;
          printTrack(m_track);
        }
        if (m_track.hits().size() > 3) {
          m_track.tagUsedHits();
          break;
        }
      }
    }
  }

}

//=========================================================================
//  Convert the local tracks to LHCb tracks
//=========================================================================
void PrPixelTracking::makeLHCbTracks() {

#ifdef DEBUG_HISTO
  unsigned int nFwd = 0; 
  unsigned int nBwd = 0;
#endif

  LHCb::Tracks* outputTracks = new LHCb::Tracks();
  put(outputTracks, m_outputLocation);
  unsigned int key = 0;
  PrPixelTracks::iterator itt;
  for (itt = m_tracks.begin(); m_tracks.end() != itt; ++itt) {
    // Skip 3-hit tracks with double-used hits.
    if ((*itt).hits().size() == 3 && (*itt).nbUnused() != 3) continue; 
    // Create a new LHCb track.
    LHCb::Track* newTrack = new LHCb::Track(key++);
    newTrack->setType(LHCb::Track::Velo);
    newTrack->setHistory(LHCb::Track::PatFastVelo);
    newTrack->setPatRecStatus(LHCb::Track::PatRecIDs);
    if (m_debug) {
      info() << "=== Store track Nb " << outputTracks->size() << endmsg;
      printTrack(*itt);
    }

    // Loop over the hits, add their LHCbIDs to the LHCb track and 
    // find the highest Z.
    double zMax = -1.e9;
    PrPixelHits::iterator ith;
    for (ith = (*itt).hits().begin(); (*itt).hits().end() != ith; ++ith) {
      newTrack->addToLhcbIDs((*ith)->id());
      if ((*ith)->z() > zMax) zMax = (*ith)->z(); 
    }
    // Decide if this is a forward or backward track.
    // Calculate Z where the track passes closest to the beam.
    const double zBeam = (*itt).zBeam();
    // Define backward as z closest to beam downstream of hits.
    const bool backward = zBeam > zMax;
    newTrack->setFlag(LHCb::Track::Backward, backward);

    // Get the state at zBeam from the straight line fit.
    LHCb::State state;
    state.setLocation(LHCb::State::ClosestToBeam);
    state.setState((*itt).state(zBeam));
    state.setCovariance((*itt).covariance(zBeam));

    // Parameters for kalmanfit scattering. calibrated on MC, shamelessly hardcoded:
    const double tx = state.tx(); const double ty = state.ty();
    const double scat2 = 1e-8 + 7e-6*(tx*tx+ty*ty) ;

    // The logic is a bit messy in the following, so I hope we got all cases right
    if (m_stateClosestToBeamKalmanFit || m_addStateFirstLastMeasurementKalmanFit) {
      // Run a K-filter with scattering to improve IP resolution
      LHCb::State upstreamstate;
      (*itt).fitKalman(upstreamstate, backward ? 1 : -1 , scat2);
      // Add this state as state at first measurement if requested
      if (m_addStateFirstLastMeasurementKalmanFit) {
        upstreamstate.setLocation(LHCb::State::FirstMeasurement);
        newTrack->addToStates(upstreamstate);
      }
      // Transport the state to the closestToBeam position
      if (m_stateClosestToBeamKalmanFit) {
        upstreamstate.setLocation(LHCb::State::ClosestToBeam);
        upstreamstate.linearTransportTo(zBeam);
        newTrack->addToStates(upstreamstate);
      }
    }
    if (!m_stateClosestToBeamKalmanFit) {
      newTrack->addToStates(state);
    }
    
    // Set state at last measurement, if requested
    if ((!backward && m_stateEndVeloKalmanFit) || m_addStateFirstLastMeasurementKalmanFit) {
      LHCb::State downstreamstate;
      (*itt).fitKalman(downstreamstate, backward ? -1 : +1 , scat2);
      if(m_addStateFirstLastMeasurementKalmanFit) {
        downstreamstate.setLocation(LHCb::State::LastMeasurement);
        newTrack->addToStates(downstreamstate);
      }
      if (m_stateEndVeloKalmanFit) {
        state = downstreamstate;
      }
    } 
    
    // Add state at end of velo
    if (!backward) {
      state.setLocation(LHCb::State::EndVelo) ;
      state.linearTransportTo(StateParameters::ZEndVelo);
      newTrack->addToStates(state);
    }

    // Set the chi2/dof
    newTrack->setNDoF(2 * ((*itt).hits().size() - 2)); 
    newTrack->setChi2PerDoF((*itt).chi2());
    // Add the LHCb track to the list.
    outputTracks->insert(newTrack);

#ifdef DEBUG_HISTO
    const unsigned int nHitsPerTrack = (*itt).hits().size();
    if (backward) {
      plot(nHitsPerTrack, "Bwd_HitsPerTrack", "Number of hits per backward track",
           0.5, 40.5, 40);
      plot(newTrack->chi2PerDoF(), "Bwd_Chi2PerTrack", "Chi2/DoF of backward tracks",
           0.0, 10.0, 50);
      plot(newTrack->pseudoRapidity(), "Bwd_EtaOfTracks", "pseudoRapidity of backward tracks",
           1.0, 6.0, 50);
      plot(newTrack->phi()*(180.0/M_PI), "Bwd_PhiOfTracks", "Phi-angle of backward tracks",
           -180.0, 180.0, 60);
      plot2D(newTrack->pseudoRapidity(), nHitsPerTrack, "Bwd_HitsPerTrackVsEta", "hits/track vs pseudoRapidity of backward tracks",
             1.0, 6.0, 0.5, 15.5, 50, 15);
      plot2D(newTrack->pseudoRapidity(), newTrack->chi2PerDoF(), "Bwd_Chi2VsEta", "Chi2/DoF vs pseudoRapidity of backward tracks",
             1.0, 6.0, 0.0, 10.0, 50, 20);
      plot2D(nHitsPerTrack, newTrack->chi2PerDoF(), "Bwd_Chi2VsHitsPerTrack", "Chi2/DoF vs hits/backward track",
             0.5, 15.5, 0.0, 10.0, 15, 20);
      nBwd++;
    } else {
      plot(nHitsPerTrack, "Fwd_HitsPerTrack", "Number of hits per forward track",
           0.5, 40.5, 40);
      plot(newTrack->chi2PerDoF(), "Fwd_Chi2PerTrack", "Chi2/DoF of forward tracks",
           0.0, 10.0, 50);
      plot(newTrack->pseudoRapidity(), "Fwd_EtaOfTracks", "pseudoRapidity of forward tracks",
           1.0, 6.0, 50);
      plot(newTrack->phi()*(180.0/M_PI), "Fwd_PhiOfTracks", "Phi-angle of forward tracks",
           -180.0, 180.0, 60);
      plot2D(newTrack->pseudoRapidity(), nHitsPerTrack, "Fwd_HitsPerTrackVsEta", "hits/track vs pseudoRapidity of forward tracks",
             1.0, 6.0, 0.5, 15.5, 50, 15);
      plot2D(newTrack->pseudoRapidity(), newTrack->chi2PerDoF(), "Fwd_Chi2VsEta", "Chi2/DoF vs pseudoRapidity of forward tracks",
             1.0, 6.0, 0.0, 10.0, 50, 20);
      plot2D(nHitsPerTrack, newTrack->chi2PerDoF(), "Fwd_Chi2VsHitsPerTrack", "Chi2/DoF vs hits/forward track",
             0.5, 15.5, 0.0, 10.0, 15, 20);
      nFwd++;
    }
#endif
  }

#ifdef DEBUG_HISTO
  plot(nFwd, "Fwd_TracksPerEvent", "Number of forward tracks per event", 0.0, 400.0, 40);
  plot(nBwd, "Bwd_TracksPerEvent", "Number of backward tracks per event", 0.0, 400.0, 40);
#endif

  m_tracks.clear();

}

//=========================================================================
//  Add hits from the specified module to the track
//=========================================================================
PrPixelHit* PrPixelTracking::bestHit(PrPixelModule* module, double xTol, double maxScatter,
                                     const PrPixelHit* h1, const PrPixelHit* h2) {
  if (module->empty()) return NULL;
  const double x1 = h1->x();
  const double y1 = h1->y();
  const double z1 = h1->z();
  const double x2 = h2->x();
  const double y2 = h2->y();
  const double z2 = h2->z();
  const double tx = (x2 - x1) / (z2 - z1);
  const double ty = (y2 - y1) / (z2 - z1); 
  // Extrapolate to the z-position of the module
  const double xGuess = x1 + tx * (module->z() - z1) - xTol;

  // If the first hit is already below this limit we can stop here.
  if (module->lastHitX() < xGuess) return NULL;
  // Do a binary search through the hits.
  unsigned int hit_start(0);
  unsigned int step(module->hits().size());
  const unsigned int module_nhits(step);
  const PrPixelHits& module_hits(module->hits());
  while (2 < step) { // quick skip of hits that are above the X-limit
    step /= 2;
    if ((module_hits[hit_start+step])->x() < xGuess) hit_start += step;
  }

  // Find the hit that matches best.
  unsigned int nFound = 0;
  double bestScatter = maxScatter;
  PrPixelHit* bestHit = NULL;
  PrPixelHit* hit(NULL);
  for (unsigned int i=hit_start; i<module_nhits; ++i) {
    hit = module_hits[i];
    const double hit_x = hit->x();
    const double hit_y = hit->y();
    const double hit_z = hit->z();
    const double dz = hit_z - z1; 
    const double xPred = x1 + tx * dz;
    const double yPred = y1 + ty * dz;
#ifdef DEBUG_HISTO
    plot((hit->x() - xPred) / xTol, "HitExtraErrPerTol", "Hit X extrapolation error / tolerance", -4.0, +4.0, 400);
#endif
    // If x-position is above prediction + tolerance, keep looking.
    if (hit_x + xTol < xPred) continue;
    // If x-position is below prediction - tolerance, stop the search.
    if (hit_x - xTol > xPred) break;
    const double dx = xPred - hit_x;
    const double dy = yPred - hit_y;
    // Skip hits outside the y-position tolerance.
    if (fabs(dy) > xTol) continue;
    const double scatter = sqrt(dx * dx + dy * dy) / fabs(hit_z - z2);
    if (scatter < bestScatter) { 
      bestHit = hit; 
      bestScatter = scatter;
    }
    if (scatter < maxScatter) ++nFound;
#ifdef DEBUG_HISTO
    plot(scatter, "HitScatter", "hit scatter [rad]", 0.0, 0.5, 500);
    plot2D(dx, dy, "Hit_dXdY", "Difference between hit and prediction in x and y [mm]", -1, 1, -1, 1, 500,500);
#endif
  }
#ifdef DEBUG_HISTO
  plot(nFound, "HitExtraCount", "Number of hits within the extrapolation window with chi2 within limits", 0.0, 10.0, 10);
#endif
  if (bestHit) {
#ifdef DEBUG_HISTO
    plot(bestScatter, "HitBestScatter", "best hit scatter [rad]",
         0.0, 0.1, 100);
#endif
    if (m_debug) printHitOnTrack(bestHit, false);
  }
  return bestHit;

}

//=========================================================================
// Debug the content of a hit
//=========================================================================
void PrPixelTracking::printHit(const PrPixelHit* hit, std::string title) {
  info() << title;
  info() << format(" module%3d x%8.3f y%8.3f z%8.2f used%2d",
                   hit->module(), hit->x(), hit->y(), hit->z(), hit->isUsed());
  if (m_debugTool) {
    LHCb::LHCbID id = hit->id();
    info() << " MC: ";
    m_debugTool->printKey(info(), id);
    if (matchKey(hit)) info() << " ***";
  }
  info() << endmsg;
}

//=========================================================================
// Print a track, with distance and chi2
//=========================================================================
void PrPixelTracking::printTrack(PrPixelTrack& track) {
  PrPixelHits::const_iterator ith;
  for (ith = track.hits().begin(); track.hits().end() != ith; ++ith) {
    printHit(*ith);
  }
}

//=========================================================================
// Print a hit on a track, with its distance.
//=========================================================================
void PrPixelTracking::printHitOnTrack(PrPixelHit* hit, bool ifMatch) {
  bool isMatching = matchKey(hit);
  isMatching = (isMatching && ifMatch) || (!isMatching && !ifMatch);
  if (isMatching) printHit(hit, "   ");
}
