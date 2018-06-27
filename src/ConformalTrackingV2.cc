#include "ConformalTrackingV2.h"

#include "ParameterParser.h"
#include "Parameters.h"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace AIDA;

// Static instantiation of the processor
ConformalTrackingV2 aConformalTrackingV2;

/*
 
 Pattern recognition code for the CLIC detector, using conformal mapping and cellular automaton
 
 */

ConformalTrackingV2::ConformalTrackingV2() : ConformalTracking("ConformalTrackingV2") {
  // Processor description
  _description = "ConformalTrackingV2 constructs tracks using a combined conformal mapping and cellular automaton approach.";

  // Input collections - tracker hits
  std::vector<std::string> inputTrackerHitCollections = {"VXDTrackerHits", "VXDEndcapTrackerHits", "ITrackerHits",
                                                         "OTrackerHits",   "ITrackerEndcapHits",   "OTrackerEndcapHits"};
  registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames", "Name of the TrackerHit input collections",
                           m_inputTrackerHitCollections, inputTrackerHitCollections);

  // Debugging collections - MC particles and relation collections
  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollectionName", "Name of the MCParticle input collection",
                          m_inputParticleCollection, std::string("MCParticle"));
  std::vector<std::string> inputRelationCollections = {"VXDTrackerHitRelations",          "VXDEndcapTrackerHitRelations",
                                                       "InnerTrackerBarrelHitsRelations", "OuterTrackerBarrelHitsRelations",
                                                       "InnerTrackerEndcapHitsRelations", "OuterTrackerEndcapHitsRelations"};
  registerInputCollections(LCIO::LCRELATION, "RelationsNames", "Name of the TrackerHit relation collections",
                           m_inputRelationCollections, inputRelationCollections);

  // Output collections - tracks
  registerOutputCollection(LCIO::TRACK, "SiTrackCollectionName", "Silicon track Collection Name", m_outputTrackCollection,
                           std::string("CATracks"));
  registerOutputCollection(LCIO::TRACKERHITPLANE, "DebugHits", "DebugHits", m_outputDebugHits, std::string("DebugHits"));

  // Parameters for tracking
  registerProcessorParameter("DebugPlots", "Plots for debugging the tracking", m_debugPlots, bool(false));
  registerProcessorParameter("RetryTooManyTracks", "retry with tightened parameters, when too many tracks are being created",
                             m_retryTooManyTracks, m_retryTooManyTracks);
  registerProcessorParameter("SortTreeResults", "sort results from kdtree search", m_sortTreeResults, m_sortTreeResults);
  registerProcessorParameter("Steps", "Steps", m_rawSteps, m_rawSteps);

  registerProcessorParameter("trackPurity", "Purity value used for checking if tracks are real or not", m_purity,
                             double(0.75));
  registerProcessorParameter("ThetaRange", "Angular range for initial cell seeding", m_thetaRange, double(0.1));
}

void ConformalTrackingV2::init() {
  ConformalTracking::init();

  parseStepParameters();
}

void ConformalTrackingV2::parseStepParameters() {
  ParameterParser::parseParameters(_stepParameters, m_rawSteps, m_inputTrackerHitCollections);
}  //parseStepParameters
