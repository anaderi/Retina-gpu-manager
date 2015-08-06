// Gaudi
#include "GaudiKernel/AlgFactory.h"
// LHCb
#include "Event/Track.h"
#include "Event/StateParameters.h"
// Local
#include "EventSerializer.h"
#include "PrPixelTracking.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrPixelTracking
//
// 2011-12-16 : Olivier Callot
//-----------------------------------------------------------------------------

DECLARE_ALGORITHM_FACTORY(PrPixelTracking)

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrPixelTracking::PrPixelTracking(const std::string &name, ISvcLocator *pSvcLocator)
    :
#ifdef DEBUG_HISTO
      GaudiTupleAlg(name, pSvcLocator),
#else
      GaudiAlgorithm(name, pSvcLocator),
#endif
      m_debugTool(NULL) {

  declareProperty("OutputTracksName",
                  m_outputLocation = LHCb::TrackLocation::Velo);
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

  // Parameters for debugging
  declareProperty("DebugToolName", m_debugToolName = "");
  declareProperty("WantedKey", m_wantedKey = -100);
  declareProperty("TimingMeasurement", m_doTiming = false);

  // Parameters for Kalman fit
  declareProperty("ClosestToBeamStateKalmanFit",
                  m_stateClosestToBeamKalmanFit = true);
  declareProperty("EndVeloStateKalmanFit", m_stateEndVeloKalmanFit = false);
  declareProperty("AddFirstLastMeasurementStatesKalmanFit",
                  m_addStateFirstLastMeasurementKalmanFit = false);

  // Parameters for 3D hit building
  declareProperty("RunOnRawBanks", m_runOnRawBanks = true);
  declareProperty("MaxClusterSize", m_maxClusterSize = PrPixel::SENSOR_PIXELS);
  declareProperty("Trigger", m_trigger = false);
}

//=============================================================================
// Destructor
//=============================================================================
PrPixelTracking::~PrPixelTracking() {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrPixelTracking::initialize() {

  info() << "--- Initializing GpuService" << endmsg;
  gpuService = svc<IGpuService>("GpuService", true);

  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) return sc;
  m_isDebug = msgLevel(MSG::DEBUG);

  // Set up the hit manager.
  m_hitManager = tool<PrPixelHitManager>("PrPixelHitManager");
  m_hitManager->setMaxClusterSize(m_maxClusterSize);
  m_hitManager->setTrigger(m_trigger);

  // Set up the debug tool.
  if ("" != m_debugToolName) m_debugTool = tool<IPrDebugTool>(m_debugToolName);
  // Set up the timing measurement.
  if (m_doTiming) {
    m_timerTool = tool<ISequencerTimerTool>("SequencerTimerTool/Timer", this);
    m_timeTotal = m_timerTool->addTimer("Total");
    m_timerTool->increaseIndent();
    m_timePrepare = m_timerTool->addTimer("Prepare");
    m_timePairs = m_timerTool->addTimer("Find by pairs");
    m_timeFinal = m_timerTool->addTimer("Store tracks");
    m_timerTool->decreaseIndent();
  }

  // use the square of the scatter to avoid calls to sqrt()
  m_maxScatter *= m_maxScatter;

#ifdef DEBUG_HISTO
  setHistoTopDir("VP/");
#endif
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
  // Build hits from clusters and sort them by global x within each module.
  if (m_doTiming) {
    m_timerTool->start(m_timeTotal);
    m_timerTool->start(m_timePrepare);
  }

  // Clean and build hits from clusters
  // From what PrPixel has:
  //  PrPixel format -> SoA (done alongside buildHits) -> AoS
  if (m_clearHits) m_hitManager->clearHits();

  if (m_runOnRawBanks) {
    m_hitManager->buildHitsFromRawBank();
  } else {
    m_hitManager->buildHits();
  }
  m_hitManager->sortByX();

  // Convert event into a format understandable by GPU
  EventSerializer      serializer;
  std::vector<uint8_t> serializedEvent;
  try {
    serializedEvent = serializer.serializeEvent(m_hitManager->hitPool(), m_hitManager->nbHits());
  } catch (const std::exception & e) {
      error() << "serialization failed; " << e.what() << endmsg;
  }
  if (serializedEvent.empty())
    info() << "--- Serialized event is empty! This should not happen!" << endmsg;

  debug() << "--- Submitting data to gpuService" << endmsg;
  debug() << "--- serializedEvent: 0x" << std::hex << (long long int)&serializedEvent[0] <<
    std::dec << ", size " << serializedEvent.size() << endmsg;

  if (m_doTiming) m_timerTool->stop(m_timePrepare);

  if (m_doTiming) m_timerTool->start(m_timePairs);

  // Perform search on the GPU
  std::vector<uint8_t> serializedTracks;

  try {
    gpuService->submitData(
        "PrPixelCudaHandler",
        serializedEvent.data(),
        serializedEvent.size(),
        allocTracks,
        &serializedTracks);
    try {
      debug() << "--- Deserializing tracks" << endmsg;
      m_tracks = serializer.deserializeTracks(serializedTracks);
    } catch (const std::exception & e) {
      error() << "deserialization failed; " << e.what() << endmsg;
    } catch (...) {
      error() << "deserialization failed; reason unknown" << endmsg;
    }
  } catch (const std::exception & e) {
    error() << "submission failed; " << e.what() << endmsg;
  } catch (...) {
    error() << "submission failed; reason unknown" << endmsg;
  }
  m_debug = true;
  if (m_doTiming) m_timerTool->stop(m_timePairs);

  // Convert temporary tracks to LHCb tracks.
  if (m_doTiming)
    m_timerTool->start(m_timeFinal);
  makeLHCbTracks();
  if (m_doTiming) {
    m_timerTool->stop(m_timeFinal);
    m_timerTool->stop(m_timeTotal);
  }

  return StatusCode::SUCCESS;
}

//=========================================================================
//  Convert the local tracks to LHCb tracks
//=========================================================================
void PrPixelTracking::makeLHCbTracks() {

  LHCb::Tracks *outputTracks = new LHCb::Tracks();
  put(outputTracks, m_outputLocation);
  unsigned int key = 0;
  auto endTracks = m_tracks.end();
  for (auto itt = m_tracks.begin(); itt != endTracks; ++itt) {
    // Skip 3-hit tracks with float-used hits.
    if ((*itt).hits().size() == 3) {
      if ((*itt).nbUnused() != 3) continue;
    } else {
      (*itt).fit();
    }
    // Create a new LHCb track.
    LHCb::Track *newTrack = new LHCb::Track(key++);
    newTrack->setType(LHCb::Track::Velo);
    newTrack->setHistory(LHCb::Track::PatFastVelo);
    newTrack->setPatRecStatus(LHCb::Track::PatRecIDs);
    m_debug = m_isDebug = true;
    if (m_debug) {
      info() << "fuck" << endmsg;
      info() << "=== Store track Nb " << outputTracks->size() << "\tnhits "
             << (*itt).hits().size() << endmsg;
      printTrack(*itt);
    }

    // Loop over the hits, add their LHCbIDs to the LHCb track and
    // find the highest Z.
    float zMax = -1.e9;
    auto endHits = (*itt).hits().cend();
    for (auto ith = (*itt).hits().cbegin(); ith != endHits; ++ith) {
      newTrack->addToLhcbIDs((*ith)->id());
      if ((*ith)->z() > zMax) zMax = (*ith)->z();
    }
    // Decide if this is a forward or backward track.
    // Calculate Z where the track passes closest to the beam.
    const float zBeam = (*itt).zBeam();
    // Define backward as z closest to beam downstream of hits.
    const bool backward = zBeam > zMax;
    newTrack->setFlag(LHCb::Track::Backward, backward);

    // Get the state at zBeam from the straight line fit.
    LHCb::State state;
    state.setLocation(LHCb::State::ClosestToBeam);
    state.setState((*itt).state(zBeam));
    state.setCovariance((*itt).covariance(zBeam));

    // Parameters for kalmanfit scattering. calibrated on MC, shamelessly
    // hardcoded:
    const float tx = state.tx();
    const float ty = state.ty();
    const float scat2 = 1e-8 + 7e-6 * (tx * tx + ty * ty);

    // The logic is a bit messy in the following, so I hope we got all cases
    // right
    if (m_stateClosestToBeamKalmanFit ||
        m_addStateFirstLastMeasurementKalmanFit) {
      // Run a K-filter with scattering to improve IP resolution
      LHCb::State upstreamstate;
      (*itt).fitKalman(upstreamstate, backward ? 1 : -1, scat2);
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
    if ((!backward && m_stateEndVeloKalmanFit) ||
        m_addStateFirstLastMeasurementKalmanFit) {
      LHCb::State downstreamstate;
      (*itt).fitKalman(downstreamstate, backward ? -1 : +1, scat2);
      if (m_addStateFirstLastMeasurementKalmanFit) {
        downstreamstate.setLocation(LHCb::State::LastMeasurement);
        newTrack->addToStates(downstreamstate);
      }
      if (m_stateEndVeloKalmanFit) {
        state = downstreamstate;
      }
    }

    // Add state at end of velo
    if (!backward) {
      state.setLocation(LHCb::State::EndVelo);
      state.linearTransportTo(StateParameters::ZEndVelo);
      newTrack->addToStates(state);
    }

    // Set the chi2/dof
    newTrack->setNDoF(2 * ((*itt).hits().size() - 2));
    newTrack->setChi2PerDoF((*itt).chi2());
    // Add the LHCb track to the list.
    outputTracks->insert(newTrack);
  }
  m_tracks.clear();

}

//=========================================================================
// Debug the content of a hit
//=========================================================================
void PrPixelTracking::printHit(const PrPixelHit *hit,
                               const std::string& title) const {
  info() << title;
  info() << format(" module%3d x%8.3f y%8.3f z%8.2f used%2d", hit->module(),
                   hit->x(), hit->y(), hit->z(), hit->isUsed());
  if (m_debugTool) {
    LHCb::LHCbID id = hit->id();
    info() << " MC: ";
    m_debugTool->printKey(info(), id);
    if (matchKey(hit)) info() << " ***";
  }
  info() << endmsg;
}

//=========================================================================
// Print all hits on a track.
//=========================================================================
void PrPixelTracking::printTrack(PrPixelTrack& track) const {
  auto end = track.hits().cend();
  for (auto it = track.hits().cbegin(); it != end; ++it) {
    printHit(*it);
  }
}

//=========================================================================
// Print a hit on a track.
//=========================================================================
void PrPixelTracking::printHitOnTrack(const PrPixelHit* hit,
                                      const bool ifMatch) const {
  bool isMatching = matchKey(hit);
  isMatching = (isMatching && ifMatch) || (!isMatching && !ifMatch);
  if (isMatching) printHit(hit, "   ");
}
