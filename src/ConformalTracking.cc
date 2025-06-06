#include "ConformalTracking.h"

#include "Cell.h"
#include "KDTree.h"
#include "KalmanTrack.h"

#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>
#include <MarlinTrk/IMarlinTrkSystem.h>
#include <MarlinTrk/MarlinTrkUtils.h>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>

#include <marlinutil/GeometryUtil.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Exceptions.h>
#include <marlin/Global.h>
#include <marlin/ProcessorEventSeeder.h>

#include <AIDA/IAnalysisFactory.h>

#include <TCanvas.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace AIDA;

// Static instantiation of the processor
ConformalTracking aConformalTracking;

/*

 Pattern recognition code for the CLIC detector, using conformal mapping and cellular automaton

 */

class TooManyTracksException : public marlin::SkipEventException {
public:
  TooManyTracksException(Processor* proc) : marlin::SkipEventException(proc) {}
};

ConformalTracking::ConformalTracking(std::string const& procName) : Processor(procName) {}

ConformalTracking::ConformalTracking() : Processor("ConformalTracking") {
  // Processor description
  _description = "ConformalTracking constructs tracks using a combined conformal mapping and cellular automaton approach.";

  // Input collections and parameters
  registerParameters();

  registerInputCollections(LCIO::TRACKERHITPLANE, "MainTrackerHitCollectionNames",
                           "Name of the TrackerHit input collections from the Main Tracker",
                           m_inputMainTrackerHitCollections,
                           {"ITrackerHits", "OTrackerHits", "ITrackerEndcapHits", "OTrackerEndcapHits"});

  registerInputCollections(LCIO::TRACKERHITPLANE, "VertexBarrelHitCollectionNames",
                           "Name of the TrackerHit input collections from the Vertex Barrel", m_inputVertexBarrelCollections,
                           {"VXDTrackerHits"});

  registerInputCollections(LCIO::TRACKERHITPLANE, "VertexEndcapHitCollectionNames",
                           "Name of the TrackerHit input collections from the Vertex Endcap", m_inputVertexEndcapCollections,
                           {"VXDEndcapTrackerHits"});

  // Parameters for tracking
  registerProcessorParameter("MaxCellAngle", "Cut on angle between two cells for cell to be valid", m_maxCellAngle,
                             double(0.035));
  registerProcessorParameter("MaxCellAngleRZ", "Cut on angle between two cells in RZ for cell to be valid", m_maxCellAngleRZ,
                             double(0.035));
  registerProcessorParameter("MaxDistance", "Maximum length of a cell (max. distance between two hits)", m_maxDistance,
                             double(0.015));
  registerProcessorParameter("HighPTCut", "pT threshold (in GeV) for enabling extendHighPT in extendTracks", m_highPTcut,
                             double(10.0));
  registerProcessorParameter("MaxChi2", "Maximum chi2/ndof for linear conformal tracks", m_chi2cut, double(300.));
  registerProcessorParameter("MinClustersOnTrack", "Minimum number of clusters to create a track in pattern recognition",
                             m_minClustersOnTrack, int(6));

  registerProcessorParameter("EnableTightCutsVertexCombined",
                             "Enabled tight cuts as first step of reconstruction in vertex b+e [TMP!!]", m_enableTCVC,
                             bool(true));
}

void ConformalTracking::registerParameters() {
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

  registerProcessorParameter("DebugPlots", "Plots for debugging the tracking", m_debugPlots, bool(false));
  registerProcessorParameter("DebugTiming", "Print out time profile", m_debugTime, bool(false));

  // Parameters for tracking
  registerProcessorParameter("RetryTooManyTracks", "retry with tightened parameters, when too many tracks are being created",
                             m_retryTooManyTracks, m_retryTooManyTracks);
  registerProcessorParameter("TooManyTracks", "Number of tracks that is considered too many", m_tooManyTracks,
                             m_tooManyTracks);
  registerProcessorParameter("SortTreeResults", "sort results from kdtree search", m_sortTreeResults, m_sortTreeResults);
  registerProcessorParameter("trackPurity", "Purity value used for checking if tracks are real or not", m_purity,
                             double(0.75));
  registerProcessorParameter("ThetaRange", "Angular range for initial cell seeding", m_thetaRange, double(0.1));
  registerProcessorParameter("MinClustersOnTrackAfterFit", "Final minimum number of track clusters",
                             m_minClustersOnTrackAfterFit, int(4));
  registerProcessorParameter("MaxHitInvertedFit", "Maximum number of track hits to try the inverted fit", m_maxHitsInvFit,
                             int(0));
}

void ConformalTracking::init() {
  // Print the initial parameters
  printParameters();

  parseStepParameters();

  // Reset counters
  m_runNumber   = 0;
  m_eventNumber = 0;

  // Set up the track fit factory
  trackFactory = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
  trackFactory->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
  trackFactory->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
  trackFactory->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, false);
  trackFactory->init();

  // Put default values for track fitting
  m_initialTrackError_d0    = 1.e6;
  m_initialTrackError_phi0  = 1.e2;
  m_initialTrackError_omega = 1.e-4;
  m_initialTrackError_z0    = 1.e6;
  m_initialTrackError_tanL  = 1.e2;
  m_maxChi2perHit           = 1.e2;

  // Get the magnetic field
  m_magneticField = MarlinUtil::getBzAtOrigin();  // z component at (0,0,0)

  // Seed hit for debug printouts. If not set later, isn't used
  debugSeed = nullptr;

  // Initialise histograms (if debug plotting on)
  if (m_debugPlots) {
    // Automatically add histograms to output root file
    AIDAProcessor::histogramFactory(this);

    // Histogram initailisation
    m_szDistribution  = new TH2F("m_szDistribution", "m_szDistribution", 20000, -100, 100, 200, -10, 10);
    m_uvDistribution  = new TH2F("m_uvDistribution", "m_uvDistribution", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    m_xyDistribution  = new TH2F("m_xyDistribution", "m_xyDistribution", 500, -1500, 1500, 500, -1500, 1500);
    m_xyzDistribution = new TH3F("m_xyzDistribution", "m_xyzDistribution", 50, 0, 100, 50, 0, 100, 100, 0, 25);

    // Histograms for neighbors parameters
    m_X      = new TH1F("m_X", "m_X", 500, -1500, 1500);
    m_Y      = new TH1F("m_Y", "m_Y", 500, -1500, 1500);
    m_Z      = new TH1F("m_Z", "m_Z", 500, -2500, 2500);
    m_neighX = new TH1F("m_neighX", "m_neighX", 500, -1500, 1500);
    m_neighY = new TH1F("m_neighY", "m_neighY", 500, -1500, 1500);
    m_neighZ = new TH1F("m_neighZ", "m_neighZ", 500, -2500, 2500);

    m_slopeZ            = new TH1F("m_slopeZ", "m_slopeZ", 1000, -100, 100);
    m_slopeZ_true       = new TH1F("m_slopeZ_true", "m_slopeZ_true", 1000, -100, 100);
    m_slopeZ_true_first = new TH1F("m_slopeZ_true_first", "m_slopeZ_true_first", 1000, -100, 100);
    m_slopeZ_vs_pt_true = new TH2F("m_diffZ_pt_true_first", "m_diffZ_pt_true_first", 2000, -500, 500, 400, 0, 100);

    // Histograms for tuning parameters (cell angle cut, cell length cut)
    m_cellAngle           = new TH1F("cellAngle", "cellAngle", 1250, 0, 0.05);
    m_cellDOCA            = new TH1F("cellDOCA", "cellDOCA", 100., 0, 0.1);
    m_cellAngleRadius     = new TH2F("cellAngleRadius", "cellAngleRadius", 400, 0, 0.04, 1000, 0, 0.04);
    m_cellLengthRadius    = new TH2F("cellLengthRadius", "cellLengthRadius", 300, 0, 0.03, 1000, 0, 0.04);
    m_cellAngleLength     = new TH2F("cellAngleLength", "cellAngleLength", 400, 0, 0.04, 300, 0, 0.03);
    m_conformalChi2       = new TH1F("conformalChi2", "conformalChi2", 100, 0, 100);
    m_conformalChi2real   = new TH1F("conformalChi2real", "conformalChi2real", 1000, 0, 1000);
    m_conformalChi2fake   = new TH1F("conformalChi2fake", "conformalChi2fake", 1000, 0, 1000);
    m_conformalChi2Purity = new TH2F("conformalChi2Purity", "conformalChi2Purity", 150, 0, 1.5, 1000, 0, 1000);

    m_conformalChi2MC        = new TH1F("conformalChi2MC", "conformalChi2MC", 1000, 0, 1000);
    m_conformalChi2PtMC      = new TH2F("conformalChi2PtMC", "conformalChi2PtMC", 1000, 0, 1000, 1000, 0, 100);
    m_conformalChi2VertexRMC = new TH2F("conformalChi2VertexRMC", "conformalChi2VertexRMC", 1000, 0, 1000, 100, 0, 100);

    m_conformalChi2SzMC   = new TH1F("conformalChi2SzMC", "conformalChi2SzMC", 1000, 0, 1000);
    m_conformalChi2SzPtMC = new TH2F("conformalChi2SzPtMC", "conformalChi2SzPtMC", 1000, 0, 1000, 1000, 0, 100);
    m_conformalChi2SzVertexRMC =
        new TH2F("conformalChi2SzVertexRMC", "conformalChi2SzVertexRMC", 1000, 0, 1000, 100, 0, 100);

    m_cellAngleMC        = new TH1F("cellAngleMC", "cellAngleMC", 1250, 0, 0.05);
    m_cellDOCAMC         = new TH1F("cellDOCAMC", "cellDOCAMC", 100., 0, 0.1);
    m_cellAngleRadiusMC  = new TH2F("cellAngleRadiusMC", "cellAngleRadiusMC", 400, 0, 0.04, 1000, 0, 0.04);
    m_cellLengthRadiusMC = new TH2F("cellLengthRadiusMC", "cellLengthRadiusMC", 300, 0, 0.03, 1000, 0, 0.04);
    m_cellAngleLengthMC  = new TH2F("cellAngleLengthMC", "cellAngleLengthMC", 400, 0, 0.04, 300, 0, 0.03);

    m_cellAngleRZMC = new TH1F("cellAngleRZMC", "cellAngleRZMC", 1250, 0, 0.05);

    // Histograms for "event display"
    m_conformalEvents       = new TH2F("conformalEvents", "conformalEvents", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    m_nonconformalEvents    = new TH2F("nonconformalEvents", "nonconformalEvents", 500, -1500, 1500, 500, -1500, 1500);
    m_conformalEventsRTheta = new TH2F("conformalEventsRTheta", "conformalEventsRTheta", 200, 0, 0.05, 632, -0.02, 6.30);
    m_conformalEventsMC     = new TH2F("conformalEventsMC", "conformalEventsMC", 1000, -0.05, 0.05, 1000, -0.05, 0.05);

    m_canvConformalEventDisplay = new TCanvas("canvConformalEventDisplay", "canvConformalEventDisplay");
    m_canvConformalEventDisplayAllCells =
        new TCanvas("canvConformalEventDisplayAllCells", "canvConformalEventDisplayAllCells");
    m_canvConformalEventDisplayAcceptedCells =
        new TCanvas("canvConformalEventDisplayAcceptedCells", "canvConformalEventDisplayAcceptedCells");
    m_canvConformalEventDisplayMC = new TCanvas("canvConformalEventDisplayMC", "canvConformalEventDisplayMC");
    m_canvConformalEventDisplayMCunreconstructed =
        new TCanvas("canvConformalEventDisplayMCunreconstructed", "canvConformalEventDisplayMCunreconstructed");
  }

  // Register this process
  Global::EVENTSEEDER->registerProcessor(this);
}

void ConformalTracking::parseStepParameters() {
  fillCollectionIndexVectors();

  int step = 0;
  // Build tracks in the vertex barrel
  Parameters initialParameters(m_vertexBarrelHits, m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut, m_minClustersOnTrack,
                               m_maxDistance, m_slopeZRange, m_highPTcut, /*highPT*/ true, /*OnlyZS*/ false,
                               /*rSearch*/ false,
                               /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                               /*combine*/ true, /*build*/ true, /*extend*/ false, /*sort*/ false);
  // Extend through the endcap
  Parameters parameters2(m_vertexEndcapHits, m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut, m_minClustersOnTrack,
                         m_maxDistance, m_slopeZRange, m_highPTcut, /*highPT*/ true, /*OnlyZS*/ false,
                         /*rSearch*/ false, /*vtt*/ true,
                         /*kalmanFitForward*/ true, step++,
                         /*combine*/ true, /*build*/ false, /*extend*/ true, /*sort*/ false);
  // Make combined vertex tracks
  Parameters parametersTCVC(m_vertexCombinedHits, m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut, m_minClustersOnTrack,
                            m_maxDistance, m_slopeZRange, m_highPTcut, /*highPT*/ true, /*OnlyZS*/ false,
                            /*rSearch*/ false, /*vtt*/ true,
                            /*kalmanFitForward*/ true, step++,
                            /*combine*/ true, /*build*/ true, /*extend*/ false, /*sort*/ false);
  // Make leftover tracks in the vertex with lower requirements
  // 1. open the cell angles
  Parameters lowerCellAngleParameters(m_vertexCombinedHits, m_maxCellAngle * 5.0, m_maxCellAngleRZ * 5.0, m_chi2cut,
                                      m_minClustersOnTrack, m_maxDistance, m_slopeZRange, m_highPTcut,
                                      /*highPT*/ true, /*OnlyZS*/ false,
                                      /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                      /*combine*/ not m_enableTCVC, /*build*/ true, /*extend*/ false, /*sort*/ false);
  // 2. open further the cell angles and increase the chi2cut
  Parameters lowerCellAngleParameters2({}, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 20.0,
                                       m_minClustersOnTrack, m_maxDistance, m_slopeZRange, m_highPTcut,
                                       /*highPT*/ true, /*OnlyZS*/ false,
                                       /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                       /*combine*/ false, /*build*/ true, /*extend*/ false, /*sort*/ false);
  // 3. min number of hits on the track = 4

  Parameters lowNumberHitsParameters({}, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 20.0,
                                     /*m_minClustersOnTrack*/ 4, m_maxDistance, m_slopeZRange, m_highPTcut,
                                     /*highPT*/ true,
                                     /*OnlyZS*/ false,
                                     /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                     /*combine*/ false, /*build*/ true, /*extend*/ false, /*sort*/ true);
  // Extend through inner and outer trackers
  Parameters trackerParameters(m_trackerHits, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 20.0,
                               /*m_minClustersOnTrack*/ 4, m_maxDistance, m_slopeZRange, /*highPTcut*/ 1.0,
                               /*highPT*/ true,
                               /*OnlyZS*/ false,
                               /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                               /*combine*/ true, /*build*/ false, /*extend*/ true, /*sort*/ false);
  // Finally reconstruct displaced tracks
  Parameters displacedParameters(m_allHits, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 10.0,
                                 /*m_minClustersOnTrack*/ 5, 0.015, m_slopeZRange, m_highPTcut,
                                 /*highPT*/ false, /*OnlyZS*/ true,
                                 /*rSearch*/ true,
                                 /*vtt*/ false, /*kalmanFitForward*/ true, step++,
                                 /*combine*/ true, /*build*/ true, /*extend*/ false, /*sort*/ false);

  _stepParameters.push_back(initialParameters);
  _stepParameters.push_back(parameters2);
  if (m_enableTCVC) {
    _stepParameters.push_back(parametersTCVC);
  }
  _stepParameters.push_back(lowerCellAngleParameters);
  _stepParameters.push_back(lowerCellAngleParameters2);
  _stepParameters.push_back(lowNumberHitsParameters);
  _stepParameters.push_back(trackerParameters);
  _stepParameters.push_back(displacedParameters);
}

// The main code, run over each event
void ConformalTracking::processEvent(LCEvent* evt) {
  //------------------------------------------------------------------------------------------------------------------
  // This pattern recognition algorithm is based on two concepts: conformal mapping and cellular automaton. Broadly
  // speaking, the 2D xy projection of all hits is transformed such that circles (helix projections) become straight
  // lines. The tracking is then considered as a 2D straight line search, using the z information to reduce
  // combinatorics.
  //
  // The hits from each input collection are transformed into conformal hit positions, with some binary trees created
  // to allow fast nearest neighbour calculations. All hits are then considered as seeds (starting from outer radius)
  // and an attempt to make cells leading back to this seed is carried out. If a long enough chain of cells can be
  // produced, this is defined as a track and the hits contained are all removed from further consideration. Once all
  // hits in a collection have been considered, the next collection of hits is added to the unused hits and the search
  // for new tracks begin again (first an attempt to extend existing tracks is performed, followed by a new search
  // using all unused hits as seeding points).
  //
  // Where several paths are possible back to the seed position, the candidate with lowest chi2/ndof is chosen.
  //------------------------------------------------------------------------------------------------------------------

  streamlog_out(DEBUG9) << "Event number: " << m_eventNumber << std::endl;

  // Set up ID decoder
  UTIL::BitField64 m_encoder(lcio::LCTrackerCellID::encoding_string());

  // Object to store all of the track hit collections passed to the pattern recognition
  std::vector<LCCollection*>        trackerHitCollections;
  std::vector<SLCRelationNavigator> relations;
  LCCollection*                     particleCollection;

  // Loop over each input collection and get the hits
  for (unsigned int collection = 0; collection < m_inputTrackerHitCollections.size(); collection++) {
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt);
    if (trackerHitCollection == 0)
      continue;
    streamlog_out(DEBUG9) << "Collection " << m_inputTrackerHitCollections[collection] << " contains "
                          << trackerHitCollection->getNumberOfElements() << " hits" << std::endl;
    trackerHitCollections.push_back(trackerHitCollection);

    // If debugging, get the relations between tracker hits and MC particle
    if (m_debugPlots) {
      // Get the collection of tracker hit relations
      LCCollection* trackerHitRelationCollection = 0;
      getCollection(trackerHitRelationCollection, m_inputRelationCollections[collection], evt);
      if (trackerHitRelationCollection == 0)
        continue;
      // Create the relations navigator
      auto relation = std::make_shared<LCRelationNavigator>(trackerHitRelationCollection);
      relations.push_back(std::move(relation));
    }
  }

  // Get the MC particle collection
  getCollection(particleCollection, m_inputParticleCollection, evt);

  // Make the output track collection
  auto trackCollection    = std::unique_ptr<LCCollectionVec>(new LCCollectionVec(LCIO::TRACK));
  auto debugHitCollection = std::unique_ptr<LCCollectionVec>(new LCCollectionVec(LCIO::TRACKERHITPLANE));
  debugHitCollection->setSubset(true);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  //------------------------------------------------------------------------------
  // Make the collection of conformal hits that will be used, with a link back to
  // the corresponding tracker hit.
  //------------------------------------------------------------------------------

  // Collections to be stored throughout the tracking
  std::map<int, SharedKDClusters>        collectionClusters;  // Conformal hits
  std::map<SKDCluster, TrackerHitPlane*> kdClusterMap;        // Their link to "real" hits
  std::map<TrackerHitPlane*, SKDCluster> conformalHits;       // The reverse link

  // Debug collections (not filled if debug off)
  std::map<MCParticle*, SharedKDClusters> particleHits;   // List of conformal hits on each MC particle
  std::map<MCParticle*, bool>             reconstructed;  // Check for MC particles
  SharedKDClusters                        debugHits;      // Debug hits for plotting
  if (m_debugPlots) {
    m_debugger.clear();
    m_debugger.setRelations(relations);
  }

  // Create the conformal hit collections for each tracker hit collection (and save the link)
  for (unsigned int collection = 0; collection < trackerHitCollections.size(); collection++) {
    // Loop over tracker hits and make conformal hit collection
    SharedKDClusters tempClusters;
    int              nHits = trackerHitCollections[collection]->getNumberOfElements();
    for (int itHit = 0; itHit < nHits; itHit++) {
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(trackerHitCollections[collection]->getElementAt(itHit));

      // Get subdetector information and check if the hit is in the barrel or endcaps
      const int celId = hit->getCellID0();
      m_encoder.setValue(celId);
      int  subdet   = m_encoder[lcio::LCTrackerCellID::subdet()];
      int  side     = m_encoder[lcio::LCTrackerCellID::side()];
      int  layer    = m_encoder[lcio::LCTrackerCellID::layer()];
      int  module   = m_encoder[lcio::LCTrackerCellID::module()];
      int  sensor   = m_encoder[lcio::LCTrackerCellID::sensor()];
      bool isEndcap = false;
      bool forward  = false;
      if (side != ILDDetID::barrel) {
        isEndcap = true;
        if (side == ILDDetID::fwd)
          forward = true;
      }

      // Make a new kd cluster
      auto kdhit = std::make_shared<KDCluster>(hit, isEndcap, forward);

      // Set the subdetector information
      kdhit->setDetectorInfo(subdet, side, layer, module, sensor);

      // Store the link between the two
      kdClusterMap[kdhit] = hit;
      conformalHits[hit]  = kdhit;
      tempClusters.push_back(kdhit);

      // Store the MC link if in debug mode
      if (m_debugPlots) {
        // Get the related simulated hit(s)
        const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects(hit);
        // Take the first hit only (TODO: this should be changed? Loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit?)
        SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
        // Get the particle belonging to that hit
        MCParticle* particle = simHit->getMCParticle();
        // Store the information (not for secondaries)
        kdParticles[kdhit] = particle;
        kdSimHits[kdhit]   = simHit;
        if (!simHit->isProducedBySecondary()) {
          particleHits[particle].push_back(kdhit);
        }
        m_debugger.registerHit(collection, kdhit, hit);
        // Draw plots for event 0
        if (m_eventNumber == 0) {
          m_conformalEvents->Fill(kdhit->getU(), kdhit->getV());
          m_nonconformalEvents->Fill(hit->getPosition()[0], hit->getPosition()[1]);
          m_conformalEventsRTheta->Fill(kdhit->getR(), kdhit->getTheta());
        }
        // Set debug hit if required
        //        if (kdhit->getX() < (-537.4) && kdhit->getX() > (-537.5) && kdhit->getY() < (-144.5) && kdhit->getY() > (-144.6))
        //          debugSeed = kdhit;
      }
    }
    collectionClusters[collection] = tempClusters;
  }

  // WHAT TO DO ABOUT THIS?? POSSIBLY MOVE DEPENDING ON MC RECONSTRUCTION (and in fact, would fit better into the check reconstruction code at present)

  // Now loop over all MC particles and make the cells connecting hits
  if (m_debugPlots) {
    int nParticles = particleCollection->getNumberOfElements();
    for (int itP = 0; itP < nParticles; itP++) {
      // Get the particle
      MCParticle* mcParticle = dynamic_cast<MCParticle*>(particleCollection->getElementAt(itP));
      // Get the vector of hits from the container
      if (particleHits.count(mcParticle) == 0)
        continue;
      SharedKDClusters trackHits = m_debugger.getAssociatedHits(mcParticle);  //particleHits[mcParticle];
      // Only make tracks with n or more hits
      if (trackHits.size() < (unsigned int)m_minClustersOnTrack)
        continue;
      // Discard low momentum particles
      double particlePt = sqrt(mcParticle->getMomentum()[0] * mcParticle->getMomentum()[0] +
                               mcParticle->getMomentum()[1] * mcParticle->getMomentum()[1]);
      // Cut on stable particles
      if (mcParticle->getGeneratorStatus() != 1)
        continue;
      // Sort the hits from larger to smaller radius
      std::sort(trackHits.begin(), trackHits.end(), sort_by_radiusKD);

      // Make a track
      auto pars    = _stepParameters[0];
      auto mcTrack = std::unique_ptr<KDTrack>(new KDTrack(pars));
      // Loop over all hits for debugging
      for (auto const& cluster : trackHits) {
        // Get the conformal clusters
        mcTrack->add(cluster);
      }

      // Fit the track and plot the chi2
      mcTrack->linearRegression();
      mcTrack->linearRegressionConformal();
      m_conformalChi2MC->Fill(mcTrack->chi2ndof());
      m_conformalChi2PtMC->Fill(mcTrack->chi2ndof(), particlePt);
      m_conformalChi2SzMC->Fill(mcTrack->chi2ndofZS());
      m_conformalChi2SzPtMC->Fill(mcTrack->chi2ndofZS(), particlePt);

      double mcVertexX = mcParticle->getVertex()[0];
      double mcVertexY = mcParticle->getVertex()[1];
      double mcVertexR = sqrt(pow(mcVertexX, 2) + pow(mcVertexY, 2));
      m_conformalChi2VertexRMC->Fill(mcTrack->chi2ndof(), mcVertexR);
      m_conformalChi2SzVertexRMC->Fill(mcTrack->chi2ndofZS(), mcVertexR);

      // Now loop over the hits and make cells - filling histograms along the way
      int nHits = trackHits.size();
      for (int itHit = 0; itHit < (nHits - 2); itHit++) {
        // Get the conformal clusters
        SKDCluster cluster0 = trackHits[itHit];
        SKDCluster cluster1 = trackHits[itHit + 1];
        SKDCluster cluster2 = trackHits[itHit + 2];

        // Make the two cells connecting these three hits
        auto cell = std::make_shared<Cell>(cluster0, cluster1);
        cell->setWeight(itHit);
        auto cell1 = std::make_shared<Cell>(cluster1, cluster2);
        cell1->setWeight(itHit + 1);

        if (itHit == 0)
          m_cellDOCAMC->Fill(cell->doca());

        // Fill the debug/tuning plots
        double angleBetweenCells   = cell->getAngle(cell1);
        double angleRZBetweenCells = cell->getAngleRZ(cell1);
        double cell0Length = sqrt(pow(cluster0->getU() - cluster1->getU(), 2) + pow(cluster0->getV() - cluster1->getV(), 2));
        double cell1Length = sqrt(pow(cluster1->getU() - cluster2->getU(), 2) + pow(cluster1->getV() - cluster2->getV(), 2));

        m_cellAngleMC->Fill(angleBetweenCells);
        m_cellAngleRadiusMC->Fill(cluster2->getR(), angleBetweenCells);
        m_cellLengthRadiusMC->Fill(cluster0->getR(), cell0Length);
        m_cellAngleLengthMC->Fill(cell1Length, angleBetweenCells);
        m_cellAngleRZMC->Fill(angleRZBetweenCells);

        // Draw cells on the first event
        if (m_eventNumber == 0) {
          // Fill the event display (hit positions)
          m_conformalEventsMC->Fill(cluster0->getU(), cluster0->getV());
          m_conformalEventsMC->Fill(cluster1->getU(), cluster1->getV());
          m_conformalEventsMC->Fill(cluster2->getU(), cluster2->getV());
          // Draw the cell lines on the event display. Use the line style to show
          // if the cells would have been cut by some of the search criteria
          m_canvConformalEventDisplayMC->cd();
          if (itHit == 0) {
            drawline(cluster0, cluster1, itHit + 1);
          }
          // Draw line style differently if the cell angle was too large
          if (angleBetweenCells > (m_maxCellAngle)) {
            drawline(cluster1, cluster2, itHit + 2, 3);
          } else {
            drawline(cluster1, cluster2, itHit + 2);
          }
        }
      }
    }
  }

  // Draw the final set of conformal hits (on top of the cell lines)
  if (m_eventNumber == 0 && m_debugPlots) {
    m_canvConformalEventDisplayMC->cd();
    m_conformalEventsMC->DrawCopy("same");
    // Draw the non-MC event display
    m_canvConformalEventDisplay->cd();
    m_conformalEvents->DrawCopy("");
    m_canvConformalEventDisplayAllCells->cd();
    m_conformalEvents->DrawCopy("");
    m_canvConformalEventDisplayAcceptedCells->cd();
    m_conformalEvents->DrawCopy("");
    m_canvConformalEventDisplayMC->cd();
    m_conformalEvents->DrawCopy("");
    m_canvConformalEventDisplayMCunreconstructed->cd();
    m_conformalEvents->DrawCopy("");
  }

  // END OF "WHAT TO DO ABOUT THIS??"

  //------------------------------------------------------------------------------
  // Now the track reconstruction strategy. Perform a sequential search, with hits
  // removed from the seeding collections once tracks have been built
  //------------------------------------------------------------------------------

  // The final vector of conformal tracks
  UniqueKDTracks   conformalTracks;
  SharedKDClusters kdClusters;
  UKDTree          nearestNeighbours = nullptr;

  for (auto const& parameters : _stepParameters) {
    try {
      runStep(kdClusters, nearestNeighbours, conformalTracks, collectionClusters, parameters);
    } catch (const TooManyTracksException& e) {
      streamlog_out(ERROR) << "Too many tracks in step: " << parameters._step << " skipping further steps" << std::endl;
      break;
    }
    streamlog_out(DEBUG9) << "STEP " << parameters._step << ": nr tracks = " << conformalTracks.size() << std::endl;
    if (streamlog_level(DEBUG9)) {
      for (auto const& confTrack : conformalTracks) {
        streamlog_out(DEBUG9) << "- Track " << &confTrack << " has " << confTrack->m_clusters.size() << " hits" << std::endl;
        for (unsigned int ht = 0; ht < confTrack->m_clusters.size(); ht++) {
          SKDCluster const& kdhit = confTrack->m_clusters.at(ht);
          streamlog_out(DEBUG9) << "-- Hit " << ht << ": [x,y,z] = [" << kdhit->getX() << ", " << kdhit->getY() << ", "
                                << kdhit->getZ() << "]" << std::endl;
        }
      }
    }
  }

  // Clean up
  nearestNeighbours.reset(nullptr);

  // Now in principle have all conformal tracks, but due to how the check for clones is performed (ish) there is a possibility
  // that clones/fakes are still present. Try to remove them by looking at overlapping hits. Turned off at the moment

  // if the conformalTracks objects needs to be used, this needs to be changed a lot
  UniqueKDTracks conformalTracksFinal = std::move(conformalTracks);

  // Sort the tracks by length, so that the longest are kept
  /*std::sort(conformalTracks.begin(),conformalTracks.end(),sort_by_length);

    for(int existingTrack=0;existingTrack<conformalTracks.size();existingTrack++){

        bool clone = false; bool saved = false;

        for(int savedTrack=0;savedTrack<conformalTracksFinal.size();savedTrack++){

        const int nOverlappingHits = overlappingHits(conformalTracks[existingTrack],conformalTracksFinal[savedTrack]);
        if( nOverlappingHits >= 2) {
            clone = true;

            // Calculate the new and existing chi2 values
            double newchi2 = (conformalTracks[existingTrack]->chi2ndofZS()*conformalTracks[existingTrack]->chi2ndofZS() + conformalTracks[existingTrack]->chi2ndof()*conformalTracks[existingTrack]->chi2ndof());
            double oldchi2 = (conformalTracksFinal[savedTrack]->chi2ndofZS()*conformalTracksFinal[savedTrack]->chi2ndofZS() + conformalTracksFinal[savedTrack]->chi2ndof()*conformalTracksFinal[savedTrack]->chi2ndof());

            double deltachi2ZS = (conformalTracks[existingTrack]->chi2ndofZS()-conformalTracksFinal[savedTrack]->chi2ndofZS());
            double deltachi2 = (conformalTracks[existingTrack]->chi2ndof()-conformalTracksFinal[savedTrack]->chi2ndof());

            // If the new track is a subtrack of an existing track, don't consider it further (already try removing bad hits from tracks
            if(nOverlappingHits == conformalTracks[existingTrack]->m_clusters.size()) break;

            // Otherwise take the longest if the delta chi2 is not too much
            else if(conformalTracks[existingTrack]->m_clusters.size() >= conformalTracksFinal[savedTrack]->m_clusters.size()){ // New track longer/equal in length

                // Increase in chi2 is too much (double)
                if( (newchi2 - oldchi2) > oldchi2) break;

                // Otherwise take it
                delete conformalTracksFinal[savedTrack];
                conformalTracksFinal[savedTrack] = conformalTracks[existingTrack]; saved = true;
            }
            else if(conformalTracks[existingTrack]->m_clusters.size() < conformalTracksFinal[savedTrack]->m_clusters.size()){ // Old track longer

                // Must improve chi2 by factor two
                if( (newchi2 - 0.5*oldchi2) > 0.) break;

                // Otherwise take it
                delete conformalTracksFinal[savedTrack];
                conformalTracksFinal[savedTrack] = conformalTracks[existingTrack]; saved = true;
            }

            break;
        }
        }

        if(!clone){
            conformalTracksFinal.push_back(conformalTracks[existingTrack]);
        }else{
            if(!saved) delete conformalTracks[existingTrack];
        }

    }
	//*/

  // Now make "real" tracks from all of the conformal tracks
  streamlog_out(DEBUG9) << "*** CA has made " << conformalTracksFinal.size()
                        << (conformalTracksFinal.size() == 1 ? " track ***" : " tracks ***") << std::endl;

  // Loop over all track candidates
  for (auto& conformalTrack : conformalTracksFinal) {
    streamlog_out(DEBUG9) << "- Fitting track " << &conformalTrack << std::endl;

    // Make the LCIO track hit vector
    EVENT::TrackerHitVec trackHits;
    for (auto const& cluster : conformalTrack->m_clusters) {
      trackHits.push_back(kdClusterMap[cluster]);
    }
    // Add kalman filtered hits
    if (conformalTrack->kalmanTrack() != nullptr) {
      KalmanTrack* kalmanTrack = conformalTrack->kalmanTrack();
      for (size_t i = 1; i < kalmanTrack->m_kalmanClusters.size(); i++) {
        SKDCluster cluster = kalmanTrack->m_kalmanClusters[i];
        trackHits.push_back(kdClusterMap[cluster]);
      }
    }

    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(), trackHits.end(), sort_by_radius);

    // Now we can make the track object and relations object, and fit the track
    auto track = std::unique_ptr<TrackImpl>(new TrackImpl);

    // First, for some reason there are 2 track objects, one which gets saved and one which is used for fitting. Don't ask...
    shared_ptr<MarlinTrk::IMarlinTrack> marlinTrack(trackFactory->createTrack());

    // Make an initial covariance matrix with very broad default values
    EVENT::FloatVec covMatrix(15, 0);             // Size 15, filled with 0s
    covMatrix[0]  = (m_initialTrackError_d0);     //sigma_d0^2
    covMatrix[2]  = (m_initialTrackError_phi0);   //sigma_phi0^2
    covMatrix[5]  = (m_initialTrackError_omega);  //sigma_omega^2
    covMatrix[9]  = (m_initialTrackError_z0);     //sigma_z0^2
    covMatrix[14] = (m_initialTrackError_tanL);   //sigma_tanl^2

    streamlog_out(DEBUG9) << " Track hits before fit = " << trackHits.size() << std::endl;

    // Try to fit
    int fitError =
        MarlinTrk::createFinalisedLCIOTrack(marlinTrack.get(), trackHits, track.get(), conformalTrack->m_kalmanFitForward,
                                            covMatrix, m_magneticField, m_maxChi2perHit);

    streamlog_out(DEBUG9) << " Fit direction " << ((conformalTrack->m_kalmanFitForward) ? "forward" : "backward")
                          << std::endl;
    streamlog_out(DEBUG9) << " Track hits after fit = " << track->getTrackerHits().size() << std::endl;

    // If the track is too short (usually less than 7 hits correspond to a vertex track)
    // the fit is tried using the inverted direction
    if (int(track->getTrackerHits().size()) < m_maxHitsInvFit || fitError != MarlinTrk::IMarlinTrack::success) {
      shared_ptr<MarlinTrk::IMarlinTrack> marlinTrack_inv(trackFactory->createTrack());
      auto                                track_inv = std::unique_ptr<TrackImpl>(new TrackImpl);

      // Try to fit on the other way
      int fitError_inv = MarlinTrk::createFinalisedLCIOTrack(marlinTrack_inv.get(), trackHits, track_inv.get(),
                                                             !conformalTrack->m_kalmanFitForward, covMatrix, m_magneticField,
                                                             m_maxChi2perHit);

      streamlog_out(DEBUG9) << " Fit direction " << ((!conformalTrack->m_kalmanFitForward) ? "forward" : "backward")
                            << std::endl;
      streamlog_out(DEBUG9) << " Track hits after inverse fit = " << track_inv->getTrackerHits().size() << std::endl;

      if (track_inv->getTrackerHits().size() > track->getTrackerHits().size()) {
        streamlog_out(DEBUG9) << " Track is replaced. " << std::endl;
        track.swap(track_inv);
        marlinTrack.swap(marlinTrack_inv);
        fitError = fitError_inv;
      } else {
        streamlog_out(DEBUG9) << " Track is not replaced. " << std::endl;
      }
    }

    // Check track quality - if fit fails chi2 will be 0. For the moment add hits by hand to any track that fails the track fit, and store it as if it were ok...
    if (fitError != MarlinTrk::IMarlinTrack::success) {
      //streamlog_out(DEBUG7) << "- Fit failed. Track has " << track->getTrackerHits().size() << " hits" << std::endl;
      streamlog_out(DEBUG9) << "- Fit fail error " << fitError << std::endl;
      continue;
    }  //*/

    // Check if track has minimum number of hits
    if (int(track->getTrackerHits().size()) < m_minClustersOnTrackAfterFit) {
      streamlog_out(DEBUG9) << "- Track has " << track->getTrackerHits().size() << " hits. The minimum required is "
                            << m_minClustersOnTrackAfterFit << std::endl;
      continue;
    }

    /*for (unsigned int p = 0; p < trackHits.size(); p++) {
      track->addHit(trackHits[p]);
    }  //*/

    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2] = trackHits.size();

    // calculate purities and check if track has been reconstructed
    if (m_debugPlots) {
      m_conformalChi2->Fill(conformalTrack->chi2ndof());
      streamlog_out(DEBUG9) << "-------------------- New TRACK --------------------" << std::endl;
      //streamlog_out(DEBUG7) << " LCIO track fit chi2 is "<<track->getChi2()<<std::endl;
      double purity = checkReal(conformalTrack, reconstructed, particleHits);
      if (purity >= m_purity) {
        m_conformalChi2real->Fill(conformalTrack->chi2ndof());
      }
      if (purity < m_purity) {
        m_conformalChi2fake->Fill(conformalTrack->chi2ndof());
      }
      m_conformalChi2Purity->Fill(purity, conformalTrack->chi2ndof());
    }

    // Push back to the output container
    trackCollection->addElement(track.release());
  }

  // Draw the cells for all produced tracks
  if (m_debugPlots && m_eventNumber == 0) {
    m_canvConformalEventDisplay->cd();
    for (auto& debugTrack : conformalTracksFinal) {
      SharedKDClusters clusters = debugTrack->m_clusters;
      std::sort(clusters.begin(), clusters.end(), sort_by_lower_radiusKD);
      for (size_t itCluster = 1; itCluster < clusters.size(); itCluster++)
        drawline(clusters[itCluster - 1], clusters[itCluster], clusters.size() - itCluster);
    }
  }

  // Draw the conformal event display hits for debugging
  if (m_debugPlots && m_eventNumber == 0) {
    m_canvConformalEventDisplay->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayAllCells->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayAcceptedCells->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayMCunreconstructed->cd();
    m_conformalEvents->DrawCopy("same");
  }
  if (m_debugPlots) {
    int nReconstructed(0), nUnreconstructed(0);
    // Additionally draw all tracks that were not reconstructed
    m_canvConformalEventDisplayMCunreconstructed->cd();
    int nParticles = particleCollection->getNumberOfElements();

    // Make the nearest neighbour tree to debug particle reconstruction issues
    SharedKDClusters kdClusters_debug;
    for (size_t i = 0; i < trackerHitCollections.size(); i++)
      kdClusters_debug.insert(kdClusters_debug.begin(), collectionClusters[i].begin(), collectionClusters[i].end());
    auto nearestNeighbours_debug = UKDTree(new KDTree(kdClusters_debug, m_thetaRange, m_sortTreeResults));

    for (int itP = 0; itP < nParticles; itP++) {
      // Get the particle
      MCParticle* mcParticle = dynamic_cast<MCParticle*>(particleCollection->getElementAt(itP));
      // Get the conformal hits
      if (particleHits.count(mcParticle) == 0)
        continue;
      SharedKDClusters mcHits = particleHits[mcParticle];
      // Cut on the number of hits
      int uniqueHits = getUniqueHits(mcHits);
      if (uniqueHits < m_minClustersOnTrack)
        continue;
      // Check if it was stable
      if (mcParticle->getGeneratorStatus() != 1)
        continue;
      // Check if it was reconstructed
      streamlog_out(DEBUG9) << "-------------------- New PARTICLE --------------------" << std::endl;
      // List the pt
      double particlePt = sqrt(mcParticle->getMomentum()[0] * mcParticle->getMomentum()[0] +
                               mcParticle->getMomentum()[1] * mcParticle->getMomentum()[1]);
      streamlog_out(DEBUG9) << "Particle pt: " << particlePt << std::endl;
      //checkReconstructionFailure(mcParticle, particleHits, nearestNeighbours_debug, _stepParameters[0]);
      if (reconstructed.count(mcParticle)) {
        nReconstructed++;
        continue;
      }
      // Draw the cells connecting the hits
      std::sort(mcHits.begin(), mcHits.end(), sort_by_radiusKD);
      for (size_t iHit = 0; iHit < (mcHits.size() - 1); iHit++) {
        drawline(mcHits[iHit], mcHits[iHit + 1], iHit + 1);
      }
      streamlog_out(DEBUG7) << "Unreconstructed particle pt: " << particlePt << std::endl;
      nUnreconstructed++;

      // Check why particles were not reconstructed
      //      checkReconstructionFailure(mcParticle, particleHits, used, nearestNeighbours);
    }
    streamlog_out(DEBUG9) << "Reconstructed " << nReconstructed << " particles out of " << nReconstructed + nUnreconstructed
                          << ". Gives efficiency "
                          << 100. * (double)nReconstructed / (double)(nReconstructed + nUnreconstructed) << "%" << std::endl;
    nearestNeighbours_debug.reset(nullptr);
  }

  // Save the output track collection
  evt->addCollection(trackCollection.release(), m_outputTrackCollection);
  evt->addCollection(debugHitCollection.release(), m_outputDebugHits);

  // Increment the event number
  m_eventNumber++;
}

void ConformalTracking::end() {
  streamlog_out(MESSAGE) << " end()  " << name() << " processed " << m_eventNumber << " events in " << m_runNumber
                         << " runs " << std::endl;

  // Write debug canvases to output file
  if (m_debugPlots) {
    m_canvConformalEventDisplay->Write();
    m_canvConformalEventDisplayAllCells->Write();
    m_canvConformalEventDisplayAcceptedCells->Write();
    m_canvConformalEventDisplayMC->Write();
    m_canvConformalEventDisplayMCunreconstructed->Write();
  }

  kdParticles.clear();
  kdSimHits.clear();
  //FIXME trackFactory is leaking Memory, but probably a MarlinTRK issue
}

//===================================
// Tracking strategies
//===================================

// Combine collections
void ConformalTracking::combineCollections(SharedKDClusters& kdClusters, UKDTree& nearestNeighbours,
                                           std::vector<int> const&                combination,
                                           std::map<int, SharedKDClusters> const& collectionClusters) {
  // Clear the input objects
  kdClusters.clear();
  nearestNeighbours.reset(nullptr);

  // Loop over all given collections
  for (unsigned int i = 0; i < combination.size(); i++) {
    // Copy the clusters to the output vector
    const SharedKDClusters& clusters = collectionClusters.at(combination[i]);
    kdClusters.insert(kdClusters.end(), clusters.begin(), clusters.end());
  }

  streamlog_out(DEBUG9) << "*** combineCollections: Collection has " << kdClusters.size() << " hits" << std::endl;

  // Sort the KDClusters from larger to smaller radius
  std::sort(kdClusters.begin(), kdClusters.end(), sort_by_radiusKD);

  // Make the binary search tree. This tree class contains two binary trees - one sorted by u-v and the other by theta
  nearestNeighbours = UKDTree(new KDTree(kdClusters, m_thetaRange, m_sortTreeResults));
}

// Take a collection of hits and try to produce tracks out of them
void ConformalTracking::buildNewTracks(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                       UKDTree& nearestNeighbours, Parameters const& parameters, bool radialSearch,
                                       bool vertexToTracker) {
  streamlog_out(DEBUG9) << "*** buildNewTracks" << std::endl;

  // Sort the input collection by radius - higher to lower if starting with the vertex detector (high R in conformal space)
  std::sort(collection.begin(), collection.end(), (vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD));

  // Loop over all hits, using each as a seed to produce a new track
  unsigned int nKDHits = collection.size();
  for (unsigned int nKDHit = 0; nKDHit < nKDHits; nKDHit++) {
    auto stopwatch_hit       = TStopwatch();
    auto stopwatch_hit_total = TStopwatch();

    // Get the kdHit and check if it has already been used (assigned to a track)
    SKDCluster kdhit = collection[nKDHit];

    streamlog_out(DEBUG9) << "Seed hit " << nKDHit << ": [x,y,z] = [" << kdhit->getX() << ", " << kdhit->getY() << ", "
                          << kdhit->getZ() << "]" << std::endl;
    if (m_debugPlots) {
      m_X->Fill(kdhit->getX());
      m_Y->Fill(kdhit->getY());
      m_Z->Fill(kdhit->getZ());
    }

    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "Starting to seed with debug cluster" << std::endl;
    if (kdhit->used()) {
      streamlog_out(DEBUG9) << "hit already used" << std::endl;
      continue;
    }

    // Debug: Plot residuals between hit and associated SimTrackerHit
    streamlog_out(DEBUG7) << "SimHit : [x,y,z] = [" << kdSimHits[kdhit]->getPosition()[0] << ", "
                          << kdSimHits[kdhit]->getPosition()[1] << ", " << kdSimHits[kdhit]->getPosition()[2] << "] "
                          << std::endl;
    //      if(kdhit->getR() < 0.003) break; // new cut - once we get to inner radius we will never make tracks. temp? TODO: make parameter? FCC (0.005 to 0.003)

    // The tracking differentiates between the first and all subsequent hits on a chain.
    // First, take the seed hit and look for sensible hits nearby to make an initial set
    // of cells. Once these are found, extrapolate the cells and look for additional hits
    // to produce a new cell along the chain. This is done mainly for speed: you ignore
    // all combinations which would be excluded when compared with the seed cell. In the
    // end we will try to produce a track starting from this seed, then begin again for
    // the next seed hit.

    // Get the initial seed cells for this hit, by looking at neighbouring hits in a given
    // radial distance. Check that they are at lower radius and further from the IP
    SharedKDClusters results;
    double           theta = kdhit->getTheta();
    // Filter already if the neighbour is used, is on the same detector layer,
    // or is in the opposite side of the detector and points inwards
    if (radialSearch)
      nearestNeighbours->allNeighboursInRadius(
          kdhit, parameters._maxDistance, results, [&kdhit, vertexToTracker](SKDCluster const& nhit) {
            if (nhit->used())
              return true;
            if (kdhit->sameLayer(nhit))
              return true;
            //not pointing in the same direction
            if (nhit->endcap() && kdhit->endcap() && (nhit->forward() != kdhit->forward()))
              return true;
            //radial conditions not met
            if ((vertexToTracker && nhit->getR() >= kdhit->getR()) || (!vertexToTracker && nhit->getR() <= kdhit->getR()))
              return true;
            return false;
          });
    else
      nearestNeighbours->allNeighboursInTheta(
          theta, m_thetaRange, results, [&kdhit, vertexToTracker](SKDCluster const& nhit) {
            if (nhit->used())
              return true;
            if (kdhit->sameLayer(nhit))
              return true;
            //not pointing in the same direction
            if (nhit->endcap() && kdhit->endcap() && (nhit->forward() != kdhit->forward()))
              return true;
            //radial conditions not met
            if ((vertexToTracker && nhit->getR() >= kdhit->getR()) || (!vertexToTracker && nhit->getR() <= kdhit->getR()))
              return true;

            return false;
          });

    if (m_debugTime)
      streamlog_out(DEBUG7) << "  Time report: Searching for " << results.size() << " neighbours took "
                            << stopwatch_hit.RealTime() * 1000 << std::scientific << " milli-seconds" << std::endl;
    stopwatch_hit.Start(true);

    streamlog_out(DEBUG9) << "Picked up " << results.size() << " neighbours from " << (radialSearch ? "radial" : "theta")
                          << " search" << std::endl;

    // Sort the neighbours by radius
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "- picked up " << results.size() << " neighbours from " << (radialSearch ? "radial" : "theta")
                            << " search" << std::endl;
    if (results.size() == 0)
      continue;
    std::sort(results.begin(), results.end(), (vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD));

    // Objects to hold cells
    SharedCells cells;
    bool        isFirst = true;

    // Make seed cells pointing inwards/outwards (conformal space)
    for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
      // Get the neighbouring hit
      SKDCluster const& nhit = results[neighbour];

      streamlog_out(DEBUG7) << "- Neighbour " << neighbour << ": [x,y,z] = [" << nhit->getX() << ", " << nhit->getY() << ", "
                            << nhit->getZ() << "]" << std::endl;

      if (!neighbourIsCompatible(nhit, kdhit, parameters._maxSlopeZ)) {
        continue;
      }

      // Check if the cell would be too long (hit very far away)
      double length2 = ((kdhit->getU() - nhit->getU()) * (kdhit->getU() - nhit->getU()) +
                        (kdhit->getV() - nhit->getV()) * (kdhit->getV() - nhit->getV()));
      if (length2 > parameters._maxDistance * parameters._maxDistance) {
        streamlog_out(DEBUG7) << "- cell between A ([x,y] = [" << kdhit->getX() << ", " << kdhit->getY()
                              << "]) and B ([x,y] = [" << nhit->getX() << ", " << nhit->getY() << "]) is too long"
                              << std::endl;
        continue;
      }

      if (m_debugPlots) {
        m_neighX->Fill(nhit->getX());
        m_neighY->Fill(nhit->getY());
        m_neighZ->Fill(nhit->getZ());

        double distanceX = nhit->getX() - kdhit->getX();
        double distanceY = nhit->getY() - kdhit->getY();
        double distance  = sqrt(distanceX * distanceX + distanceY * distanceY);
        double deltaZ    = nhit->getZ() - kdhit->getZ();
        double slopeZ    = deltaZ / distance;
        m_slopeZ->Fill(slopeZ);

        // Debug using the seed hit and the associated SimTrackerHit
        double simpt;
        if (kdParticles[kdhit] == kdParticles[nhit] && kdhit != nhit) {
          simpt = sqrt(kdParticles[kdhit]->getMomentum()[0] * kdParticles[kdhit]->getMomentum()[0] +
                       kdParticles[kdhit]->getMomentum()[1] * kdParticles[kdhit]->getMomentum()[1]);
          streamlog_out(DEBUG7) << "- They were produced by the same MCParticle with pt = " << simpt << std::endl;
          streamlog_out(DEBUG7) << "- SimHit : [x,y,z] = [" << kdSimHits[nhit]->getPosition()[0] << ", "
                                << kdSimHits[nhit]->getPosition()[1] << ", " << kdSimHits[nhit]->getPosition()[2] << "] "
                                << std::endl;
          streamlog_out(DEBUG7) << "- Delta : [x,y,z] = [" << nhit->getX() - kdhit->getX() << ", "
                                << nhit->getY() - kdhit->getY() << ", " << nhit->getZ() - kdhit->getZ() << "] " << std::endl;

          m_slopeZ_true->Fill(slopeZ);
          m_slopeZ_vs_pt_true->Fill(slopeZ, simpt);

          if (isFirst) {
            streamlog_out(DEBUG7) << "- is first" << std::endl;
            m_slopeZ_true_first->Fill(slopeZ);
            isFirst = false;
          }
        }
        //Debugging: uncomment in the case you want to create seeds only with hits belonging to the same MCParticle
        // else {
        //  continue;
        //}
      }

      // Create the new seed cell
      cells.emplace_back(std::make_shared<Cell>(kdhit, nhit));
      auto const& cell = cells.back();

      //      if (cell->doca() > 0.01) {
      //        continue;
      //      }

      streamlog_out(DEBUG7) << "- made cell between A ([x,y] = [" << kdhit->getX() << ", " << kdhit->getY()
                            << "]) and B ([x,y] = [" << nhit->getX() << ", " << nhit->getY() << "])" << std::endl;
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG7) << "- made cell with neighbour " << neighbour << " at " << nhit->getU() << "," << nhit->getV()
                              << std::endl;

      // Debug plotting
      if (m_debugPlots) {
        m_cellDOCA->Fill(cell->doca());
        if (m_eventNumber == 0) {
          m_canvConformalEventDisplayAllCells->cd();
          drawline(kdhit, nhit, 1);
        }
      }
    }

    if (m_debugTime)
      streamlog_out(DEBUG7) << "  Time report: Making " << cells.size() << " seed cells took "
                            << stopwatch_hit.RealTime() * 1000 << std::scientific << " milli-seconds" << std::endl;
    stopwatch_hit.Start(true);

    streamlog_out(DEBUG9) << "Produced " << cells.size() << " seed cells from seed hit A ([x,y] = [" << kdhit->getX() << ", "
                          << kdhit->getY() << "])" << std::endl;
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "- produced " << cells.size() << " seed cells" << std::endl;

    // No seed cells produced
    if (cells.size() == 0)
      continue;

    // All seed cells have been created, now try create all "downstream" cells until no more can be added
    SharedKDClusters debugHits;
    if (debugSeed && kdhit == debugSeed) {
      extendSeedCells(cells, nearestNeighbours, true, debugHits, parameters, vertexToTracker);
    } else {
      extendSeedCells(cells, nearestNeighbours, true, debugHits, parameters, vertexToTracker);
    }

    if (m_debugTime)
      streamlog_out(DEBUG7) << "  Time report: Extending " << cells.size() << " seed cells took "
                            << stopwatch_hit.RealTime() * 1000 << std::scientific << " milli-seconds" << std::endl;
    stopwatch_hit.Start(true);

    streamlog_out(DEBUG9) << "After extension, have " << cells.size() << " cells from seed hit A ([x,y] = [" << kdhit->getX()
                          << ", " << kdhit->getY() << "])" << std::endl;
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "- after extension, have " << cells.size() << " cells" << std::endl;

    // Now have all cells stemming from this seed hit. If it is possible to produce a track (ie. cells with depth X) then we will now...
    //      if(depth < (m_minClustersOnTrack-1)) continue; // TODO: check if this is correct

    // We create all acceptable tracks by looping over all cells with high enough weight to create
    // a track and trace their route back to the seed hit. We then have to choose the best candidate
    // at the end (by minimum chi2 of a linear fit)
    std::map<SCell, bool> usedCells;
    std::map<SCell, bool> usedCells2;
    UniqueKDTracks        cellTracks;

    // Sort Cells from highest to lowest weight
    std::sort(cells.begin(), cells.end(), sort_by_cellWeight);

    // Create tracks by following a path along cells
    streamlog_out(DEBUG9) << "Create tracks by following path along cells. Loop over cells, sorted by weight" << std::endl;

    int nCells = cells.size();
    for (int itCell = 0; itCell < nCells; itCell++) {
      streamlog_out(DEBUG7) << "- Cell " << itCell << " between A ([x,y] = [" << cells[itCell]->getStart()->getX() << ", "
                            << cells[itCell]->getStart()->getY() << "]) and B ([x,y] = [" << cells[itCell]->getEnd()->getX()
                            << ", " << cells[itCell]->getEnd()->getY() << "]) has weight " << cells[itCell]->getWeight()
                            << std::endl;
      // Check if this cell has already been used
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG7) << "-- looking at cell " << itCell << std::endl;
      if (usedCells.count(cells[itCell])) {
        streamlog_out(DEBUG7) << "-- used cell" << std::endl;
        continue;
      }
      // Check if this cell could produce a track (is on a long enough chain)

      if (cells[itCell]->getWeight() < (parameters._minClustersOnTrack - 2)) {
        streamlog_out(DEBUG7) << "-- cell can not produce a track: weight < (minClustersOnTrack - 2)" << std::endl;
        break;
      }
      // Produce all tracks leading back to the seed hit from this cell
      UniqueCellularTracks candidateTracks;
      UniqueCellularTracks candidateTracksTemp;
      createTracksNew(candidateTracksTemp, cells[itCell],
                      usedCells2);  // Move back to using used cells here? With low chi2/ndof?

      copy_if(std::make_move_iterator(candidateTracksTemp.begin()), std::make_move_iterator(candidateTracksTemp.end()),
              std::back_inserter(candidateTracks), [parameters](UcellularTrack const& track) {
                return (int(track->size()) >= (parameters._minClustersOnTrack - 1));
              });
      candidateTracksTemp.clear();

      streamlog_out(DEBUG7) << "- From cell, produced " << candidateTracks.size() << " candidate tracks" << std::endl;
      if (streamlog_level(DEBUG9)) {
        for (auto& candidateTrack : candidateTracks) {
          streamlog_out(DEBUG7) << "--  track is made of " << candidateTrack->size() << " cells " << std::endl;
          for (unsigned int candidateCell = 0; candidateCell < candidateTrack->size(); candidateCell++) {
            streamlog_out(DEBUG7) << "--- cell between A ([x,y] = [" << (*candidateTrack)[candidateCell]->getStart()->getX()
                                  << ", " << (*candidateTrack)[candidateCell]->getStart()->getY() << "]) and B ([x,y] = ["
                                  << (*candidateTrack)[candidateCell]->getEnd()->getX() << ", "
                                  << (*candidateTrack)[candidateCell]->getEnd()->getY() << "])" << std::endl;
          }
        }
      }
      // Debug plotting
      if (m_debugPlots && m_eventNumber == 2) {
        m_canvConformalEventDisplayAcceptedCells->cd();
        for (auto& track : candidateTracks) {
          for (unsigned int trackCell = 0; trackCell < track->size(); trackCell++) {
            drawline((*track)[trackCell]->getStart(), (*track)[trackCell]->getEnd(), track->size() - trackCell);
          }
        }
      }

      // Look at the candidate tracks and fit them (+pick the best chi2/ndof)
      //       streamlog_out(DEBUG7)<<"- produced "<<candidateTracks.size()<<" candidate tracks"<<std::endl;
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG7) << "- produced " << candidateTracks.size() << " candidate tracks" << std::endl;
      if (candidateTracks.size() == 0)
        continue;
      std::vector<double> chi2ndof;
      UniqueKDTracks      bestTracks;

      // Temporary check of how many track candidates should not strictly have been allowed
      //checkUnallowedTracks(candidateTracks, parameters);

      getFittedTracks(bestTracks, candidateTracks, usedCells,
                      parameters);  // Returns all tracks at the moment, not lowest chi2 CHANGE ME

      // Store track(s) for later
      cellTracks.insert(cellTracks.end(), std::make_move_iterator(bestTracks.begin()),
                        std::make_move_iterator(bestTracks.end()));
    }

    streamlog_out(DEBUG9) << "Final number of fitted tracks to this seed hit: " << cellTracks.size() << std::endl;
    if (m_debugTime)
      streamlog_out(DEBUG7) << "  Time report: " << cellTracks.size() << " tracks reconstructed from cells took "
                            << stopwatch_hit.RealTime() * 1000 << std::scientific << " milli-seconds" << std::endl;
    stopwatch_hit.Start(true);

    // All tracks leading back to the seed hit have now been found. Decide which are the feasible candidates (may be more than 1)
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "== final number of candidate tracks to this seed hit: " << cellTracks.size() << std::endl;
    if (cellTracks.size() == 0) {
      continue;
    }

    UniqueKDTracks bestTracks = std::move(cellTracks);  //CHANGE ME - temp to give all tracks
    //      UniqueKDTracks bestTracks = getLowestChi2(cellTracks,cellTracksChi2ndof,chi2ndof);
    //     streamlog_out(DEBUG7)<<"== final number of stored tracks to this seed hit: "<<bestTracks.size()<<std::endl;
    if (debugSeed && kdhit == debugSeed) {
      streamlog_out(DEBUG7) << "== final number of stored tracks to this seed hit: " << bestTracks.size() << std::endl;
      for (size_t itBest = 0; itBest < bestTracks.size(); itBest++)
        streamlog_out(DEBUG7) << "- track " << itBest << " has chi2/ndof " << bestTracks[itBest]->chi2ndof() << std::endl;
    }

    // Could now think to do the full helix fit and apply a chi2 cut. TODO

    // Store the final CA tracks. First sort them by length, so that if we produce a clone with just 1 hit missing, the larger track is taken

    streamlog_out(DEBUG9) << "Now cut on chi2 and treate the clones. Loop over tracks, sorted by length" << std::endl;

    std::sort(bestTracks.begin(), bestTracks.end(), sort_by_length);
    for (auto& bestTrack : bestTracks) {
      bool bestTrackUsed = false;

      // Cut on chi2
      //        if(collection != trackerHitCollections.size() && (bestTracks[itTrack]->chi2ndof() > parameters._chi2cut || bestTracks[itTrack]->chi2ndofZS() > parameters._chi2cut)) {
      double chi2cut = parameters._chi2cut;
      //if (bestTracks[itTrack]->pt() < 5.)
      //  chi2cut = 1000.;

      if ((parameters._onlyZSchi2cut && bestTrack->chi2ndofZS() > chi2cut) ||
          (!parameters._onlyZSchi2cut && (bestTrack->chi2ndof() > chi2cut || bestTrack->chi2ndofZS() > chi2cut))) {
        streamlog_out(DEBUG9) << "- Track has chi2 too large" << std::endl;
        bestTrack.reset();
        continue;
      }

      // temp cut to attack fake tracks
      //if( parameters._onlyZSchi2cut && bestTracks[itTrack]->clusters().size() == 4 && bestTracks[itTrack]->chi2ndof() > 100. ){
      //  bestTracks[itTrack].reset();
      //  continue;
      //}

      // Check if the new track is a clone
      bool clone = false;

      for (auto& conformalTrack : conformalTracks) {
        const unsigned int nOverlappingHits = overlappingHits(bestTrack, conformalTrack);
        if (nOverlappingHits >= 2) {
          clone = true;
          streamlog_out(DEBUG9) << "- Track is a clone" << std::endl;

          // Calculate the new and existing chi2 values
          double newchi2 = bestTrack->chi2ndofZS() + bestTrack->chi2ndof();
          double oldchi2 = conformalTrack->chi2ndofZS() + conformalTrack->chi2ndof();

          // If the new track is an existing track + segment, take the new track
          if (nOverlappingHits == conformalTrack->m_clusters.size()) {
            conformalTrack = std::move(bestTrack);
            bestTrackUsed  = true;
            streamlog_out(DEBUG9) << "- New = existing + segment. Replaced existing with new" << std::endl;
          }

          // If the new track is a subtrack of an existing track, don't consider it further (already try removing bad hits from tracks

          else if (nOverlappingHits == bestTrack->m_clusters.size()) {
            streamlog_out(DEBUG9) << "- New + segment = existing. New ignored" << std::endl;
            break;
          }
          // Otherwise take the longest
          else if (bestTrack->m_clusters.size() == conformalTrack->m_clusters.size()) {  // New track equal in length
            if (newchi2 > oldchi2) {
              streamlog_out(DEBUG9) << "- New equal length. Worse chi2:" << newchi2 << std::endl;
              break;
            }
            // Take it
            conformalTrack = std::move(bestTrack);
            bestTrackUsed  = true;
            streamlog_out(DEBUG9) << "- New equal. Better chi2 Replaced existing with new" << std::endl;

          } else if (bestTrack->m_clusters.size() > conformalTrack->m_clusters.size()) {  // New track longer

            // Take it
            conformalTrack = std::move(bestTrack);
            bestTrackUsed  = true;
            streamlog_out(DEBUG9) << "- New longer. Replaced existing with new" << std::endl;

          } else if (bestTrack->m_clusters.size() < conformalTrack->m_clusters.size()) {  // Old track longer
            streamlog_out(DEBUG9) << "- Old longer. New ignored" << std::endl;
            break;
          }
          break;
        }
      }

      // If not a clone, save the new track
      if (!clone) {
        bestTrackUsed = true;

        streamlog_out(DEBUG9) << "- Track is not a clone. Pushing back best track with chi2/ndof " << bestTrack->chi2ndof()
                              << " and hits: " << std::endl;

        if (debugSeed && kdhit == debugSeed) {
          streamlog_out(DEBUG7) << "== Pushing back best track with chi2/ndof " << bestTrack->chi2ndof() << std::endl;
        } else {
          streamlog_out(DEBUG7) << "Pushing back best track with chi2/ndof " << bestTrack->chi2ndof() << std::endl;
        }

        for (unsigned int cluster = 0; cluster < bestTrack->m_clusters.size() && streamlog_level(DEBUG9); cluster++) {
          streamlog_out(DEBUG9) << "-- Hit " << cluster << ": [x,y] = [" << bestTrack->m_clusters.at(cluster)->getX() << ", "
                                << bestTrack->m_clusters.at(cluster)->getY() << "]" << std::endl;
        }

        conformalTracks.push_back(std::move(bestTrack));
      }

      if (not bestTrackUsed) {
        bestTrack.reset();
      }

    }  // end for besttracks
    if (m_debugTime)
      streamlog_out(DEBUG7) << "  Time report: Sort best tracks took " << stopwatch_hit.RealTime() * 1000 << std::scientific
                            << " milli-seconds" << std::endl;
    if (m_debugTime)
      streamlog_out(DEBUG7) << " Time report: Total time for seed hit " << nKDHit << " = "
                            << stopwatch_hit_total.RealTime() * 1000 << std::scientific << " milli-seconds" << std::endl;
  }
}

// Check that the neighbour has the z slope in the required range compared to the seed
bool ConformalTracking::neighbourIsCompatible(const SKDCluster& neighbourHit, const SKDCluster& seedHit,
                                              const double slopeZRange) {
  const double distanceX(neighbourHit->getX() - seedHit->getX());
  const double distanceY(neighbourHit->getY() - seedHit->getY());
  const double distance(sqrt(distanceX * distanceX + distanceY * distanceY));
  const double deltaZ(neighbourHit->getZ() - seedHit->getZ());
  if (fabs(deltaZ / distance) > slopeZRange) {
    streamlog_out(DEBUG7) << "- z condition not met" << std::endl;
    if (debugSeed && seedHit == debugSeed)
      streamlog_out(DEBUG7) << "- z condition not met" << std::endl;
    return false;
  }

  return true;
}
// Take a collection of tracks and try to extend them into the collection of clusters passed.
void ConformalTracking::extendTracks(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                     UKDTree& nearestNeighbours, Parameters const& parameters) {
  // Loop over all current tracks. At the moment this is a "stupid" algorithm: it will simply try to add every
  // hit in the collection to every track, and keep the ones thta have a good chi2. In fact, it will extrapolate
  // the track and do a nearest neighbours search, but this seemed to fail for some reason, TODO!

  int nTracks = conformalTracks.size();
  streamlog_out(DEBUG7) << "EXTENDING " << nTracks << " tracks, into " << collection.size() << " hits" << std::endl;
  if (collection.size() == 0)
    return;

  // Sort the hit collection by layer
  std::sort(collection.begin(), collection.end(), sort_by_layer);

  // First extend high pt tracks
  extendTracksPerLayer(conformalTracks, collection, nearestNeighbours, parameters, parameters._vertexToTracker);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }

  //index just for debug
  int debug_idxTrack = 0;

  for (auto& track : conformalTracks) {
    // Make sure that track hits are ordered from largest to smallest radius
    std::sort(track->m_clusters.begin(), track->m_clusters.end(), sort_by_radiusKD);

    // TODO: are kalman tracks necessary?
    // Make sure that all tracks have a kalman track attached to them
    //    if(track->kalmanTrack() == nullptr){
    //      KalmanTrack* kalmanTrack = new KalmanTrack(track);
    //      track->setKalmanTrack(kalmanTrack);
    //    }

    // Get the associated MC particle
    MCParticle* associatedParticle = nullptr;
    if (m_debugPlots) {
      associatedParticle = m_debugger.getAssociatedParticle(track);
      if (associatedParticle == nullptr)
        streamlog_out(DEBUG7) << "- nullptr particle!" << std::endl;

      // Check the track pt estimate
      TLorentzVector mc_helper;
      mc_helper.SetPxPyPzE(associatedParticle->getMomentum()[0], associatedParticle->getMomentum()[1],
                           associatedParticle->getMomentum()[2], associatedParticle->getEnergy());

      streamlog_out(DEBUG7) << "- extending track " << debug_idxTrack << " with pt = " << mc_helper.Pt()
                            << ". pt estimate: " << track->pt() << " chi2/ndof " << track->chi2ndof() << " and chi2/ndof ZS "
                            << track->chi2ndofZS() << std::endl;
    }
    debug_idxTrack++;

    // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
    int  nclusters = track->m_clusters.size();
    auto seedCell  = std::make_shared<Cell>(track->m_clusters[nclusters - 2], track->m_clusters[nclusters - 1]);

    // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
    //double searchDistance = parameters._maxDistance;
    //KDCluster* fakeHit = extrapolateCell(seedCell, searchDistance / 2.);  // TODO: make this search a function of radius
    //SharedKDClusters results2;
    //nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25 * searchDistance / 2., results2); //TEMP while not using
    //std::sort(results2.begin(), results2.end(), sort_by_radiusKD);
    //delete fakeHit;

    // Loop over all hits found and check if any have a sufficiently good delta Chi2
    SharedKDClusters goodHits;
    SKDCluster       bestCluster = nullptr;
    double           bestChi2    = 0.;

    unsigned int nKDHits = collection.size();
    //streamlog_out(DEBUG7)<<"STARTING"<<std::endl;
    //      for(int newHit=0;newHit<results2.size();newHit++){
    for (unsigned int nKDHit = 0; nKDHit < nKDHits; nKDHit++) {
      // Get the kdHit and check if it has already been used (assigned to a track)
      SKDCluster kdhit      = collection[nKDHit];
      bool       associated = false;
      if (m_debugPlots) {
        associated = m_debugger.isAssociated(kdhit, associatedParticle);
        if (associated)
          streamlog_out(DEBUG7) << "-- hit " << nKDHit << " belongs to this track" << std::endl;
      }
      //streamlog_out(DEBUG7)<<"Detector "<<kdhit->getSubdetector()<<", layer "<<kdhit->getLayer()<<", side "<<kdhit->getSide()<<std::endl;

      // If this hit is on a new sensor, then add the hit from the previous sensor and start anew
      if (bestCluster != nullptr && !(kdhit->sameSensor(bestCluster))) {
        bestCluster->used(true);
        track->add(bestCluster);
        track->linearRegression();
        track->linearRegressionConformal();
        bestCluster = nullptr;
      }

      // Don't reuse hits
      if (kdhit->used()) {
        if (associated)
          streamlog_out(DEBUG7) << "used" << std::endl;
        continue;
      }

      // Don't pick up hits in the opposite side of the detector
      if ((track->m_clusters[nclusters - 1]->getZ() > 0. && kdhit->getZ() < 0.) ||
          (track->m_clusters[nclusters - 1]->getZ() < 0. && kdhit->getZ() > 0.)) {
        if (associated)
          streamlog_out(DEBUG7) << "opposite side of detector" << std::endl;
        continue;
      }

      // First check that the hit is not wildly away from the track (make cell and check angle)
      // SCell extensionCell = std::make_shared<Cell>(track->m_clusters[0],results2[newHit]);
      Cell   extensionCell(track->m_clusters[nclusters - 1], kdhit);
      double cellAngle   = seedCell->getAngle(extensionCell);
      double cellAngleRZ = seedCell->getAngleRZ(extensionCell);
      if (cellAngle > 3. * parameters._maxCellAngle || cellAngleRZ > 3. * parameters._maxCellAngleRZ) {
        if (associated)
          streamlog_out(DEBUG7) << "-- killed by cell angle cut" << std::endl;
        continue;
      }

      // Now fit the track with the new hit and check the increase in chi2
      double deltaChi2(0.), deltaChi2zs(0.);

      fitWithPoint(*(track), kdhit, deltaChi2, deltaChi2zs);  //track->deltaChi2(results2[newHit]);

      if (associated) {
        streamlog_out(DEBUG7) << "-- hit was fitted and has a delta chi2 of " << deltaChi2 << " and delta chi2zs of "
                              << deltaChi2zs << std::endl;
      }

      // We have an estimate of the pT here, could use it in the chi2 criteria
      double chi2cut = parameters._chi2cut;
      //if (track->pt() < 5.)
      //  chi2cut = 1000.;

      if (deltaChi2 > chi2cut || deltaChi2zs > chi2cut)
        continue;

      bool onSameSensor = false;
      for (auto const& clusterOnTrack : track->m_clusters) {
        if (kdhit->sameSensor(clusterOnTrack)) {
          onSameSensor = true;
          break;
        }
      }
      if (onSameSensor)
        continue;

      if (associated)
        streamlog_out(DEBUG7) << "-- valid candidate!" << std::endl;

      if (bestCluster == nullptr || deltaChi2 < bestChi2) {
        bestCluster = kdhit;
        bestChi2    = deltaChi2;
      }
    }

    if (bestCluster != nullptr) {
      bestCluster->used(true);
      track->add(bestCluster);
      track->linearRegression();
      track->linearRegressionConformal();
      bestCluster = nullptr;
    }

  }  // End of loop over tracks
}

// Extend seed cells
void ConformalTracking::extendSeedCells(SharedCells& cells, UKDTree& nearestNeighbours, bool extendingTrack,
                                        const SharedKDClusters& /*debugHits*/, Parameters const& parameters,
                                        bool vertexToTracker) {
  streamlog_out(DEBUG8) << "***** extendSeedCells" << std::endl;

  unsigned int nCells   = 0;
  int          depth    = 0;
  int          startPos = 0;

  // Keep track of existing cells in case there are branches in the track
  std::map<SKDCluster, SharedCells> existingCells;

  streamlog_out(DEBUG8) << "Extending " << cells.size()
                        << " cells. Start with seed cells (depth 0) and then move forward (depth > 0)." << std::endl;
  // Try to create all "downstream" cells until no more can be added
  while (cells.size() != nCells) {
    // Extend all cells with depth N. In the next iteration, look at cells with depth N+1
    nCells = cells.size();
    streamlog_out(DEBUG8) << "Depth = " << depth << std::endl;

    for (unsigned int itCell = startPos; itCell < nCells; itCell++) {
      streamlog_out(DEBUG8) << "- Extend cell " << itCell << ": A ([x,y] = [" << cells[itCell]->getStart()->getX() << ", "
                            << cells[itCell]->getStart()->getY() << "]) - B ([x,y] = [" << cells[itCell]->getEnd()->getX()
                            << ", " << cells[itCell]->getEnd()->getY() << "])" << std::endl;
      // Get the end point of the cell (to search for neighbouring hits to form new cells connected to this one)
      SKDCluster const& hit            = cells[itCell]->getEnd();
      double            searchDistance = parameters._maxDistance;  //hit->getR();
      if (searchDistance > hit->getR())
        searchDistance = 1.2 * hit->getR();

      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      SKDCluster const& fakeHit =
          extrapolateCell(cells[itCell], searchDistance / 2.);  // TODO: make this search a function of radius
      SharedKDClusters results;
      nearestNeighbours->allNeighboursInRadius(
          fakeHit, 0.625 * searchDistance, results, [&hit, vertexToTracker](SKDCluster const& nhit) {
            if (nhit->used())
              return true;
            if (hit->sameLayer(nhit))
              return true;
            if (nhit->endcap() && hit->endcap() && (nhit->forward() != hit->forward()))
              return true;
            if ((vertexToTracker && nhit->getR() >= hit->getR()) || (!vertexToTracker && nhit->getR() <= hit->getR()))
              return true;
            return false;
          });

      streamlog_out(DEBUG8) << "- Found " << results.size() << " neighbours from cell extrapolation " << std::endl;
      if (extendingTrack) {
        streamlog_out(DEBUG7) << "Extrapolating cell " << itCell << " from u,v: " << hit->getU() << "," << hit->getV()
                              << std::endl;
        streamlog_out(DEBUG7) << "Extrapolated hit at u,v: " << fakeHit->getU() << "," << fakeHit->getV() << std::endl;
        streamlog_out(DEBUG7) << "Found " << results.size() << " neighbours from cell extrapolation" << std::endl;
      }

      // Make new cells pointing inwards
      for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
        // Get the neighbouring hit
        SKDCluster const& nhit = results[neighbour];

        streamlog_out(DEBUG8) << "-- Neighbour " << neighbour << ": [x,y] = [" << nhit->getX() << ", " << nhit->getY() << "]"
                              << std::endl;
        if (extendingTrack)
          streamlog_out(DEBUG7) << "looking at neighbour " << neighbour << " at u,v: " << nhit->getU() << "," << nhit->getV()
                                << std::endl;

        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
        //        if(used.count(nhit)){if(extendingTrack)streamlog_out(DEBUG7)<<"- used"<<std::endl; continue;}

        if (nhit->used()) {
          streamlog_out(DEBUG8) << "-- used" << std::endl;
          continue;
        }
        // if (hit->sameLayer(nhit)) {
        //   if (extendingTrack)
        //     streamlog_out(DEBUG7) << "- same layer" << std::endl;
        //   continue;
        // }
        // if (nhit->endcap() && hit->endcap() && (nhit->forward() != hit->forward()))
        //   continue;
        // if ((vertexToTracker && nhit->getR() >= hit->getR()) || (!vertexToTracker && nhit->getR() <= hit->getR())) {
        //   if (extendingTrack)
        //     streamlog_out(DEBUG7) << "- " << (vertexToTracker ? "higher radius" : "lower radius") << std::endl;
        //   continue;
        // }

        // Check if this cell already exists (rejoining branch) FIXME - allows rejoining a branch without checking cell angles
        auto const& existingCellsForHit = existingCells.find(hit);
        if (existingCellsForHit != existingCells.end()) {
          bool alreadyExists = false;
          for (auto const& existingCell : existingCellsForHit->second) {
            if (existingCell->getEnd() == nhit) {
              streamlog_out(DEBUG8) << "-- cell A ([x,y] = [" << hit->getX() << ", " << hit->getY() << "]) - B ([x,y] = ["
                                    << nhit->getX() << ", " << nhit->getY() << "]) already exists" << std::endl;
              alreadyExists = true;

              // Check if cell angle is too large to rejoin
              if (cells[itCell]->getAngle(existingCell) > parameters._maxCellAngle ||

                  cells[itCell]->getAngleRZ(existingCell) > parameters._maxCellAngleRZ) {
                streamlog_out(DEBUG8) << "-- cell A ([x,y] = [" << hit->getX() << ", " << hit->getY() << "]) - B ([x,y] = ["
                                      << nhit->getX() << ", " << nhit->getY() << "]) angle too large" << std::endl;
                continue;
              }
              // Otherwise add the path
              cells[itCell]->setTo(existingCell);
              existingCell->setFrom(cells[itCell]);
              updateCell(existingCell);
            }
          }
          if (alreadyExists)
            continue;
        }

        // Make the new cell
        Cell cell(hit, nhit);

        streamlog_out(DEBUG8) << "-- made new cell A ([x,y] = [" << hit->getX() << ", " << hit->getY() << "]) - B ([x,y] = ["
                              << nhit->getX() << ", " << nhit->getY() << "]) " << std::endl;
        if (extendingTrack)
          streamlog_out(DEBUG7) << "- make new cell" << std::endl;

        // Check if the new cell is compatible with the previous cell (angle between the two is acceptable)
        //        if( cells[itCell]->getAngle(cell) > (parameters._maxCellAngle*exp(-0.001/nhit->getR())) ){
        if (cells[itCell]->getAngle(cell) > parameters._maxCellAngle ||
            cells[itCell]->getAngleRZ(cell) > parameters._maxCellAngleRZ) {
          // Debug plotting
          //          if(m_debugPlots && m_eventNumber == 0){
          //            m_canvConformalEventDisplayAllCells->cd();
          //            drawline(hit,nhit,cells[itCell]->getWeight()+2,3);
          //          }

          streamlog_out(DEBUG8) << "-- cell angle too large. Discarded! " << std::endl;

          if (extendingTrack)
            streamlog_out(DEBUG7) << "-- discarded!" << std::endl;

          continue;
        }

        // Set the information about which cell this new cell is attached to and store it
        cells.emplace_back(std::make_shared<Cell>(std::move(cell)));
        auto const& scell = cells.back();
        existingCells[hit].push_back(scell);
        scell->setFrom(cells[itCell]);
        cells[itCell]->setTo(scell);

        // Debug plotting
        //        if(m_debugPlots && m_eventNumber == 0){
        //          m_canvConformalEventDisplayAllCells->cd();
        //          drawline(hit,nhit,cells[itCell]->getWeight()+2);
        //        }
      }

      // Finished adding new cells to this cell
    }

    // All new cells added at this depth
    startPos = nCells;
    depth++;
  }
  streamlog_out(DEBUG8) << "extendSeedCells *****" << std::endl;

  // No more downstream cells can be added
}

void ConformalTracking::extendHighPT(UniqueKDTracks& conformalTracks, SharedKDClusters& /*collection*/,
                                     UKDTree& nearestNeighbours, Parameters const& parameters, bool /*radialSearch*/) {
  // Set the cell angle criteria for high pt tracks
  Parameters highPTParameters(parameters);
  highPTParameters._maxCellAngle   = 0.005;
  highPTParameters._maxCellAngleRZ = 0.005;

  // Loop over all tracks
  for (auto& track : conformalTracks) {
    // Only look at high pt tracks
    if (track->pt() < highPTParameters._highPTcut)
      continue;

    streamlog_out(DEBUG7) << "Trying to extend high pt track with pt estimate: " << track->pt() << std::endl;
    streamlog_out(DEBUG7) << "- track has size " << track->m_clusters.size() << std::endl;

    // Make sure that track hits are ordered from largest to smallest radius
    std::sort(track->m_clusters.begin(), track->m_clusters.end(), sort_by_radiusKD);

    // Make seed cell from first two hits
    int  nclusters = track->m_clusters.size();
    auto seedCell  = std::make_shared<Cell>(track->m_clusters[nclusters - 2], track->m_clusters[nclusters - 1]);

    SharedKDClusters const& trackHits = track->m_clusters;
    for (size_t th = 0; th < trackHits.size(); th++) {
      streamlog_out(DEBUG7) << "Hit " << th << " u = " << trackHits[th]->getU() << " v = " << trackHits[th]->getV()
                            << " x = " << trackHits[th]->getX() << " y = " << trackHits[th]->getY()
                            << " z = " << trackHits[th]->getZ() << std::endl;
    }

    // Extend downstream
    SharedCells cells;
    cells.push_back(seedCell);
    SharedKDClusters debugHits;
    extendSeedCells(cells, nearestNeighbours, false, debugHits, highPTParameters);

    // Now make all possible track extensions
    std::map<SCell, bool> usedCells2;
    UniqueCellularTracks  extensions;

    if (cells.size() == 1)
      continue;

    // Sort Cells from highest to lowest weight
    std::sort(cells.begin(), cells.end(), sort_by_cellWeight);

    // Create tracks by following a path along cells
    int nCells = cells.size();
    for (int itCell = 0; itCell < nCells; itCell++) {
      // Produce all tracks leading back to the seed hit from this cell
      UniqueCellularTracks extensionsTemp;
      createTracksNew(extensionsTemp, cells[itCell], usedCells2);  // Move back to using used cells here? With low chi2/ndof?

      // Save them
      extensions.insert(extensions.end(), std::make_move_iterator(extensionsTemp.begin()),
                        std::make_move_iterator(extensionsTemp.end()));
    }

    // Now have all possible extensions. Get the lowest chi2
    double   lowestChi2ndof = 10000000.;
    UKDTrack lowestChi2Track;
    for (auto& extension : extensions) {
      // Get the extension, add it to the original hit
      SharedKDClusters hits = track->clusters();

      if (extension->size() < 2)
        continue;
      //hits.push_back((*extension)[extension->size()-1]->getEnd());
      for (int cell = (extension->size() - 2); cell >= 0; cell--) {
        //hits.push_back((*extension)[cell]->getStart());
        hits.push_back((*extension)[cell]->getEnd());
      }

      /*if (itTrack == 0) {
        for (size_t th = 0; th < hits.size(); th++) {
          streamlog_out(DEBUG7) << "Extension 0. Hit " << th << " u = " << hits[th]->getU() << " v = " << hits[th]->getV()
                                << " x = " << hits[th]->getX() << " y = " << hits[th]->getY() << " z = " << hits[th]->getZ()
                                << std::endl;
        }
      }*/

      UKDTrack extendedTrack = std::unique_ptr<KDTrack>(new KDTrack(parameters));
      for (auto const& hit : hits)
        extendedTrack->add(hit);
      extendedTrack->linearRegression(true);
      extendedTrack->linearRegressionConformal();
      double chi2ndof = sqrt(extendedTrack->chi2ndof() * extendedTrack->chi2ndof() +
                             extendedTrack->chi2ndofZS() * extendedTrack->chi2ndofZS());
      if (chi2ndof < lowestChi2ndof) {
        lowestChi2ndof  = chi2ndof;
        lowestChi2Track = std::move(extendedTrack);
      }
    }

    // If this is good enough then add it to the track
    if (lowestChi2Track == nullptr)
      continue;

    streamlog_out(DEBUG7) << "Track chi2 is " << track->chi2ndof() << ", zs is " << track->chi2ndofZS() << std::endl;
    streamlog_out(DEBUG7) << "Best chi2 is " << lowestChi2Track->chi2ndof() << ", zs is " << lowestChi2Track->chi2ndofZS()
                          << std::endl;

    double existingChi2 = sqrt(track->chi2ndof() * track->chi2ndof() + track->chi2ndofZS() * track->chi2ndofZS());
    if (lowestChi2ndof < (existingChi2 + 10.)) {
      //delete track;
      track = std::move(lowestChi2Track);
    }
    streamlog_out(DEBUG7) << "- track finishes with size " << track->m_clusters.size() << std::endl;
  }

  return;
}

// Take a collection of tracks and find the best extension layer by layer in the passed collection of clusters
void ConformalTracking::extendTracksPerLayer(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                             UKDTree& nearestNeighbours, Parameters const& parameters,
                                             bool vertexToTracker) {
  streamlog_out(DEBUG9) << "*** extendTracksPerLayer: extending " << conformalTracks.size()
                        << " tracks layer by layer, into " << collection.size() << " hits" << std::endl;

  if (collection.size() == 0)
    return;

  int nTracks = conformalTracks.size();
  streamlog_out(DEBUG9) << "Total number of tracks = " << nTracks << std::endl;

  //index just for debug
  int debug_idxTrack = 0;

  // Sort the input collection by layer
  std::sort(collection.begin(), collection.end(), (vertexToTracker ? sort_by_layer : sort_by_lower_layer));

  // Loop over all tracks
  for (auto& track : conformalTracks) {
    streamlog_out(DEBUG9) << "Track " << debug_idxTrack << std::endl;
    debug_idxTrack++;

    // Only look at high pt tracks
    if (track->pt() < parameters._highPTcut)
      continue;

    // Make sure that the hits are ordered in KDradius
    std::sort(track->m_clusters.begin(), track->m_clusters.end(),
              (vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD));

    // Make the cell from the two hits of the track from which to extend
    int nclusters = track->m_clusters.size();
    streamlog_out(DEBUG9) << "Track has " << nclusters << " hits" << std::endl;
    for (int i = 0; i < nclusters; i++) {
      streamlog_out(DEBUG9) << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
                            << track->m_clusters.at(i)->getY() << "]" << std::endl;
    }
    auto seedCell = std::make_shared<Cell>(track->m_clusters[nclusters - 2], track->m_clusters[nclusters - 1]);
    streamlog_out(DEBUG9) << "Seed cell A ([x,y] = [" << track->m_clusters[nclusters - 2]->getX() << ","
                          << track->m_clusters[nclusters - 2]->getY() << "]) - B ([x,y] = ["
                          << track->m_clusters[nclusters - 1]->getX() << "," << track->m_clusters[nclusters - 1]->getY()
                          << "])" << std::endl;

    // Initialize variables for nearest neighbours search and subdetector layers
    SharedKDClusters results;
    bool             loop           = true;
    int              extendInSubdet = 0;
    int              extendInLayer  = 0;
    int              final_subdet   = 0;
    int              final_layer    = 0;
    bool             firstLayer     = true;
    bool             skipLayer      = false;
    SKDCluster       expectedHit    = nullptr;

    // Start the extension into input collection, layer by layer
    do {
      streamlog_out(DEBUG9) << "- Start extension" << std::endl;
      // Find nearest neighbours in theta and sort them by layer
      SKDCluster const& kdhit = seedCell->getEnd();
      streamlog_out(DEBUG9) << "- Endpoint seed cell: [x,y] = [" << kdhit->getX() << ", " << kdhit->getY() << "]; "
                            << "]; r = " << kdhit->getR() << "; [subdet,layer] = [" << kdhit->getSubdetector() << ", "
                            << kdhit->getLayer() << "]" << std::endl;
      double theta = kdhit->getTheta();
      nearestNeighbours->allNeighboursInTheta(
          theta, m_thetaRange * 4, results, [&kdhit, vertexToTracker](SKDCluster const& nhit) {
            //if same subdet, take only the hits within two layers
            if (nhit->getSubdetector() == kdhit->getSubdetector()) {
              if ((vertexToTracker && nhit->getLayer() > (kdhit->getLayer() + 2)) ||
                  (!vertexToTracker && nhit->getLayer() < (kdhit->getLayer() - 2))) {
                return true;
              }
            }

            if ((vertexToTracker && nhit->getR() > kdhit->getR()) || (!vertexToTracker && nhit->getR() < kdhit->getR()))
              return true;
            return false;
          });
      //nearestNeighbours->allNeighboursInRadius(kdhit, parameters._maxDistance, results);
      streamlog_out(DEBUG9) << "- Found " << results.size() << " neighbours. " << std::endl;
      if (streamlog_level(DEBUG9)) {
        for (auto const& neighbour : results) {
          streamlog_out(DEBUG9) << "-- Neighbour from allNeighboursInTheta : [x,y] = [" << neighbour->getX() << ", "
                                << neighbour->getY() << "]; r = " << neighbour->getR() << "; [subdet,layer] = ["
                                << neighbour->getSubdetector() << ", " << neighbour->getLayer() << "]" << std::endl;
        }
      }

      if (results.size() == 0) {
        loop = false;
        continue;
      }

      std::sort(results.begin(), results.end(), (vertexToTracker ? sort_by_layer : sort_by_lower_layer));
      // Get the final values of subdet and layer to stop the loop at the (vertexToTracker? outermost : innermost) layer with neighbours
      final_subdet = results.back()->getSubdetector();
      final_layer  = results.back()->getLayer();
      streamlog_out(DEBUG9) << "- Final layer with neighbours: [subdet,layer] = [" << final_subdet << ", " << final_layer
                            << "]" << std::endl;

      // If no hit was found on a layer, we use expectedHit as current position
      // However, the nearest neighbour search is always performed from the last real hit (kdhit)
      const SKDCluster& hitOnCurrentLayer = skipLayer ? expectedHit : kdhit;

      // Set the subdetector layer in which to extend -- for this it is important that the hits are sorted by layer before
      // First layer is based on first neighbour
      // Next layers are based on next wrt current
      for (auto const& neighbour : results) {
        streamlog_out(DEBUG9) << "-- Neighbour : [x,y] = [" << neighbour->getX() << ", " << neighbour->getY()
                              << "]; [subdet,layer] = [" << neighbour->getSubdetector() << ", " << neighbour->getLayer()
                              << "]" << std::endl;

        if (firstLayer) {  // first neighbour
          extendInSubdet = results.front()->getSubdetector();
          extendInLayer  = results.front()->getLayer();
          streamlog_out(DEBUG9) << "-- firstLayer) [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer << "]"
                                << std::endl;
          firstLayer = false;
          break;
        } else if (neighbour->getSubdetector() ==
                   hitOnCurrentLayer->getSubdetector()) {  // next layer is in same subdetector
          if ((vertexToTracker && (neighbour->getLayer() > hitOnCurrentLayer->getLayer())) ||
              (!vertexToTracker && (neighbour->getLayer() < hitOnCurrentLayer->getLayer()))) {
            extendInSubdet = neighbour->getSubdetector();
            extendInLayer  = neighbour->getLayer();
            streamlog_out(DEBUG9) << "-- sameSub,diffLayer) [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer
                                  << "]" << std::endl;
            break;
          }
        } else if ((vertexToTracker && (neighbour->getSubdetector() > hitOnCurrentLayer->getSubdetector())) ||
                   (!vertexToTracker && (neighbour->getSubdetector() <
                                         hitOnCurrentLayer->getSubdetector()))) {  // next layer is in next subdetector
          extendInSubdet = neighbour->getSubdetector();
          extendInLayer  = neighbour->getLayer();
          streamlog_out(DEBUG9) << "-- diffSub) [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer << "]"
                                << std::endl;
          break;
        } else
          continue;
      }

      // Set the condition to end the loop
      if (extendInSubdet == final_subdet && extendInLayer == final_layer) {
        streamlog_out(DEBUG9) << "- Ending the loop. Reached final subdet layer" << std::endl;
        loop = false;
      }

      // Initialize variables for choosing the best neighbour in layer
      SKDCluster              bestCluster          = nullptr;
      map<SKDCluster, double> bestClustersWithChi2 = {};
      vector<SKDCluster>      bestClusters         = {};
      double                  bestChi2             = std::numeric_limits<double>::max();
      double                  chi2                 = track->chi2ndof();
      double                  chi2zs               = track->chi2ndofZS();
      KDTrack                 tempTrack            = *track;

      // Loop over neighbours
      for (auto const& neighbour : results) {
        // only in the layer in which to extend
        if (neighbour->getSubdetector() == extendInSubdet && neighbour->getLayer() == extendInLayer) {
          // Store the hit in case no bestCluster will be found in this layer
          expectedHit = neighbour;

          streamlog_out(DEBUG9) << "-- Looking at neighbour " << neighbour << ": [x,y] = [" << neighbour->getX() << ", "
                                << neighbour->getY() << "]" << std::endl;

          // Check that the hit has not been used
          if (neighbour->used()) {
            streamlog_out(DEBUG9) << "-- used" << std::endl;
            continue;
          }

          // Check that the hit is not in the opposite side of the detector (if endcap)
          if (neighbour->endcap() && kdhit->endcap() && (neighbour->forward() != kdhit->forward())) {
            streamlog_out(DEBUG9) << "-- opposite side of detector" << std::endl;
            continue;
          }

          // Check that radial conditions are met
          if ((vertexToTracker && neighbour->getR() >= kdhit->getR()) ||
              (!vertexToTracker && neighbour->getR() <= kdhit->getR())) {
            streamlog_out(DEBUG9) << "-- radial conditions not met" << std::endl;
            continue;
          }

          // Check that the hit is not wildly away from the track (make cell and check angle)
          Cell   extensionCell(track->m_clusters[nclusters - 1], neighbour);
          double cellAngle      = seedCell->getAngle(extensionCell);
          double cellAngleRZ    = seedCell->getAngleRZ(extensionCell);
          double maxCellAngle   = parameters._maxCellAngle;
          double maxCellAngleRZ = parameters._maxCellAngleRZ;

          if (cellAngle > maxCellAngle || cellAngleRZ > maxCellAngleRZ) {
            streamlog_out(DEBUG9) << "-- killed by cell angle cut" << std::endl;
            continue;
          }

          // Now fit the track with the new hit and check the increase in chi2 - use a tempTrack object: add hit, fit, get chi2, remove hit
          double deltaChi2(0.), deltaChi2zs(0.);
          streamlog_out(DEBUG9) << "-- tempTrack has " << tempTrack.m_clusters.size() << " hits " << std::endl;
          tempTrack.add(neighbour);
          tempTrack.linearRegression();
          tempTrack.linearRegressionConformal();
          double newchi2   = tempTrack.chi2ndof();
          double newchi2zs = tempTrack.chi2ndofZS();
          streamlog_out(DEBUG9) << "-- tempTrack has now " << tempTrack.m_clusters.size() << " hits " << std::endl;
          deltaChi2   = newchi2 - chi2;
          deltaChi2zs = newchi2zs - chi2zs;
          streamlog_out(DEBUG9) << "-- hit was fitted and has a delta chi2 of " << deltaChi2 << " and delta chi2zs of "
                                << deltaChi2zs << std::endl;
          tempTrack.remove(tempTrack.m_clusters.size() - 1);
          streamlog_out(DEBUG9) << "-- tempTrack has now " << tempTrack.m_clusters.size() << " hits " << std::endl;

          double chi2cut = parameters._chi2cut;
          if (deltaChi2 > chi2cut || deltaChi2zs > chi2cut) {
            streamlog_out(DEBUG9) << "-- killed by chi2 cut" << std::endl;
            continue;
          }
          streamlog_out(DEBUG9) << "-- valid candidate" << std::endl;

          bestClustersWithChi2[neighbour] = deltaChi2;

          // bestCluster still empty - fill it with the first candidate
          // otherwise fill it with the one with best chi2
          if (bestCluster == nullptr || deltaChi2 < bestChi2) {
            bestCluster = neighbour;
            bestChi2    = deltaChi2;
          } else {
            continue;
          }

        }  // end if on the extension layer

      }  // end loop on neighbours

      streamlog_out(DEBUG9) << "-- this seed cells has " << bestClustersWithChi2.size() << " good candidates." << std::endl;

      //put the best cluster already found
      if (bestCluster != nullptr) {
        bestClusters.push_back(bestCluster);
      }

      //put all the other cluster with a similar chi2 - probably duplicates
      float chi2window = 100.0;
      if (bestClustersWithChi2.size() > 0) {
        for (auto clu : bestClustersWithChi2) {
          streamlog_out(DEBUG9) << "-- Best cluster candidate: [x,y] = [" << clu.first->getX() << ", " << clu.first->getY()
                                << "]; r = " << clu.first->getR() << "; radius = " << clu.first->getRadius() << std::endl;
          if (clu.second < bestChi2 + chi2window && clu.second > bestChi2 - chi2window) {
            if (!clu.first->sameSensor(bestCluster))
              bestClusters.push_back(clu.first);
            else
              streamlog_out(DEBUG9) << "-- current candidate and best one are in the same sensor " << std::endl;
          }
        }
      }
      streamlog_out(DEBUG9) << "-- this seed cells will be updated with " << bestClusters.size() << " candidates."
                            << std::endl;
      std::sort(bestClusters.begin(), bestClusters.end(), (vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD));

      // If bestClusters have been found in this layer, add them to the track, update the seed cell and reset
      if (!bestClusters.empty()) {
        skipLayer = false;
        nclusters = track->m_clusters.size();
        streamlog_out(DEBUG9) << "- nclusters = " << nclusters << std::endl;
        for (int i = 0; i < nclusters; i++) {
          streamlog_out(DEBUG9) << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
                                << track->m_clusters.at(i)->getY() << "]; r = " << track->m_clusters.at(i)->getR()
                                << "; radius = " << track->m_clusters.at(i)->getRadius() << std::endl;
        }
        //create new cell with last track cluster and last bestCluster
        //Best clusters are already ordered depending on vertexToTracker bool
        seedCell = std::make_shared<Cell>(track->m_clusters[nclusters - 1], bestClusters.at(bestClusters.size() - 1));

        for (auto bestClu : bestClusters) {
          streamlog_out(DEBUG9) << "- Found bestCluster [x,y] = [" << bestClu->getX() << ", " << bestClu->getY()
                                << "]; r = " << bestClu->getR() << std::endl;
          track->add(bestClu);
          track->linearRegression();
          track->linearRegressionConformal();
          bestClu->used(true);
          nclusters = track->m_clusters.size();
          streamlog_out(DEBUG9) << "- nclusters = " << nclusters << std::endl;
          for (int i = 0; i < nclusters; i++) {
            streamlog_out(DEBUG9) << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
                                  << track->m_clusters.at(i)->getY() << "]; r = " << track->m_clusters.at(i)->getR()
                                  << "; radius = " << track->m_clusters.at(i)->getRadius() << std::endl;
          }
        }
      }
      // If not bestCluster has been found in this layer, make cell with the expected hit (from extrapolation) and increment the missing hit count
      else {
        streamlog_out(DEBUG9) << "- Found no bestCluster for [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer
                              << "]" << std::endl;
        skipLayer = true;
      }
      streamlog_out(DEBUG9) << (skipLayer ? "- Still" : "- Updated") << " seed cell A ([x,y] = ["
                            << seedCell->getStart()->getX() << ", " << seedCell->getStart()->getY() << "]) - B ([x,y] = ["
                            << seedCell->getEnd()->getX() << ", " << seedCell->getEnd()->getY() << ")]" << std::endl;

      // Clear the neighbours tree
      results.clear();
      bestClustersWithChi2.clear();
      bestClusters.clear();

    } while (loop);  // end of track extension

    streamlog_out(DEBUG9) << "Track ends with " << track->m_clusters.size() << " hits" << std::endl;
    if (streamlog_level(DEBUG9)) {
      for (unsigned int i = 0; i < track->m_clusters.size(); i++) {
        streamlog_out(DEBUG9) << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
                              << track->m_clusters.at(i)->getY() << "], [subdet,layer] = ["
                              << track->m_clusters.at(i)->getSubdetector() << ", " << track->m_clusters.at(i)->getLayer()
                              << "]" << std::endl;
      }
    }

  }  // end loop on tracks
}

//===================================
// Cellular Track Functions
//===================================

// New test at creating cellular tracks. In this variant, don't worry about clones etc, give all possible routes back to the seed cell. Then cut
// on number of clusters on each track, and pass back (good tracks to then be decided based on best chi2
void ConformalTracking::createTracksNew(UniqueCellularTracks& finalcellularTracks, SCell& seedCell,
                                        std::map<SCell, bool>& /*usedCells*/) {
  streamlog_out(DEBUG8) << "***** createTracksNew" << std::endl;

  // Final container to be returned
  UniqueCellularTracks cellularTracks;

  // Make the first cellular track using the seed cell
  auto seedTrack = std::unique_ptr<cellularTrack>(new cellularTrack());
  seedTrack->push_back(seedCell);
  cellularTracks.push_back(std::move(seedTrack));

  streamlog_out(DEBUG8) << "Follow all paths from higher weighted cells back to the seed cell" << std::endl;
  // Now start to follow all paths back from this seed cell
  // While there are still tracks that are not finished (last cell weight 0), keep following their path
  while (toBeUpdated(cellularTracks)) {
    //   streamlog_out(DEBUG7)<<"== Updating "<<cellularTracks.size()<<" tracks"<<std::endl;
    // Loop over all (currently existing) tracks
    int nTracks = cellularTracks.size();

    streamlog_out(DEBUG8) << "Going to create " << nTracks << " tracks" << std::endl;
    if (nTracks > 10000) {
      streamlog_out(WARNING) << "WARNING: Going to create " << nTracks << " tracks " << std::endl;
    }
    if (nTracks > m_tooManyTracks) {
      streamlog_out(ERROR) << "Too many tracks (" << nTracks << " > " << m_tooManyTracks
                           << " are going to be created, tightening parameters" << std::endl;
      throw TooManyTracksException(this);
    }
    for (int itTrack = 0; itTrack < nTracks; itTrack++) {
      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom().size() == 0) {
        streamlog_out(DEBUG8) << "- Cellular track " << itTrack << " is finished " << std::endl;
        //       streamlog_out(DEBUG7)<<"-- Track "<<itTrack<<" is finished"<<std::endl;
        continue;
      }

      // While there is only one path leading from this cell, follow that path
      SCell cell = cellularTracks[itTrack]->back();
      //     streamlog_out(DEBUG7)<<"-- Track "<<itTrack<<" has "<<(*(cell->getFrom())).size()<<" cells attached to the end of it"<<std::endl;
      //      while(cell->getWeight() > 0 && (*(cell->getFrom())).size() == 1){
      while (cell->getFrom().size() == 1) {
        streamlog_out(DEBUG8) << "- Cellular track " << itTrack << " is a simple extension" << std::endl;
        //       streamlog_out(DEBUG7)<<"- simple extension"<<std::endl;
        // Get the cell that it attaches to
        auto parentCell = SCell(cell->getFrom()[0]);
        streamlog_out(DEBUG8) << "- Added parent cell A ([x,y] = [" << parentCell->getStart()->getX() << ", "
                              << parentCell->getStart()->getY() << "]) - B ([x,y] = [" << cell->getEnd()->getX() << ", "
                              << cell->getEnd()->getY() << "])" << std::endl;
        // Attach it to the track and continue
        cellularTracks[itTrack]->push_back(parentCell);
        cell = parentCell;
      }

      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom().size() == 0)
        continue;

      // If the weight is != 0 and there is more than one path to follow, branch the track (create a new one for each path)
      int nBranches = cell->getFrom().size();

      streamlog_out(DEBUG8) << "- Cellular track " << itTrack << " has more than one extension" << std::endl;

      //     streamlog_out(DEBUG7)<<"- making "<<nBranches<<" branches"<<std::endl;

      // For each additional branch make a new track
      for (int itBranch = 1; itBranch < nBranches; itBranch++) {
        streamlog_out(DEBUG8) << "-- Cellular track " << itTrack << ", extension " << itBranch << std::endl;
        auto branchedTrack      = std::unique_ptr<cellularTrack>(new cellularTrack(*(cellularTracks[itTrack].get())));
        auto branchedParentCell = SCell(cell->getFrom()[itBranch]);
        streamlog_out(DEBUG8) << "-- Added branched parent cell A ([x,y] = [" << branchedParentCell->getStart()->getX()
                              << ", " << branchedParentCell->getStart()->getY() << "]) - B ([x,y] = ["
                              << cell->getEnd()->getX() << ", " << cell->getEnd()->getY() << "])" << std::endl;
        branchedTrack->push_back(std::move(branchedParentCell));
        cellularTracks.push_back(std::move(branchedTrack));
      }

      // Keep the existing track for the first branch
      cellularTracks[itTrack]->push_back(SCell(cell->getFrom()[0]));
    }
  }

  finalcellularTracks.insert(finalcellularTracks.end(), std::make_move_iterator(cellularTracks.begin()),
                             std::make_move_iterator(cellularTracks.end()));

  streamlog_out(DEBUG8) << "Number of finalcellularTracks = " << finalcellularTracks.size() << std::endl;
  if (streamlog_level(DEBUG8)) {
    for (auto& cellTrack : finalcellularTracks) {
      streamlog_out(DEBUG8) << "- Finalcelltrack is made of " << cellTrack->size() << " cells " << std::endl;
      for (unsigned int finalcell = 0; finalcell < cellTrack->size(); finalcell++) {
        streamlog_out(DEBUG8) << "-- Cell A ([x,y] = [" << (*cellTrack)[finalcell]->getStart()->getX() << ", "
                              << (*cellTrack)[finalcell]->getStart()->getY() << "]) - B ([x,y] = ["
                              << (*cellTrack)[finalcell]->getEnd()->getX() << ", "
                              << (*cellTrack)[finalcell]->getEnd()->getY() << "])" << std::endl;
      }
    }
  }
  streamlog_out(DEBUG8) << "createTracksNew *****" << std::endl;

  return;
}

// Check if any of the tracks in a collection still have to be updated
bool ConformalTracking::toBeUpdated(UniqueCellularTracks const& cellularTracks) {
  for (auto const& aCellularTrack : cellularTracks)
    if (aCellularTrack->back()->getFrom().size() > 0) {
      return true;
    }
  return false;
}

// Given a list of connected cells (so-called cellular tracks), return the candidate(s) with lowest chi2/degrees of freedom.
// Attempt to remove hits to see if there is a significant improvement in the chi2/ndof, to protect against noise hits causing
// a good track to be discarded. If several candidates have low chi2/ndof and are not clones (limited sharing of hits) then
// return all of them. Given that there is no material scattering taken into account this helps retain low pt tracks, which
// may have worse chi2/ndof than ghosts/real tracks with an additional unrelated hit from the low pt track.
void ConformalTracking::getFittedTracks(UniqueKDTracks& finalTracks, UniqueCellularTracks&      candidateTracks,
                                        std::map<SCell, bool>& /*usedCells*/, Parameters const& parameters) {
  streamlog_out(DEBUG8) << "***** getFittedTracks" << std::endl;

  // Make a container for all tracks being considered, initialise variables
  UniqueKDTracks trackContainer;
  //  std::vector<double> trackChi2ndofs;

  // Loop over all candidate tracks and do an inital fit to get the track angle (needed to calculate the
  // hit errors for the error-weighted fit)
  for (auto& candidateTrack : candidateTracks) {
    // If there are not enough hits on the track, ignore it
    if (int(candidateTrack->size()) < (parameters._minClustersOnTrack - 2)) {
      streamlog_out(DEBUG8) << "- Track " << candidateTrack.get() << ": not enough hits on track." << std::endl;
      candidateTrack.reset();
      continue;
    }

    // Make the fitting object. TGraphErrors used for 2D error-weighted fitting
    UKDTrack track = std::unique_ptr<KDTrack>(new KDTrack(parameters));

    // Loop over all hits and add them to the fitter (and track)
    double     npoints = 0.;
    SKDCluster kdStart = (*candidateTrack)[0]->getEnd();
    track->add(kdStart);
    npoints++;

    for (auto& trackCell : (*candidateTrack)) {
      SKDCluster kdEnd = trackCell->getStart();
      track->add(kdEnd);
      npoints++;
    }
    track->linearRegression(parameters._highPTfit);
    track->linearRegressionConformal();

    double chi2ndofTOT = track->chi2ndof() + track->chi2ndofZS();

    // We try to see if there are spurious hits causing the chi2 to be very large. This would cause us to throw away
    // good tracks with perhaps just a single bad hit. Try to remove each hit and see if fitting without it causes a
    // significant improvement in chi2/ndof

    // Loop over each hit (starting at the back, since we will use the 'erase' function to get rid of them)
    // and see if removing it improves the chi2/ndof
    double removed = 0;
    if (chi2ndofTOT > parameters._chi2cut &&
        chi2ndofTOT <
            parameters
                ._chi2cut) {  //CHANGE ME?? Upper limit to get rid of really terrible tracks (temp lower changed from 0 to parameters._chi2cut)

      streamlog_out(DEBUG8) << "- Track " << candidateTrack.get() << ": chi2ndofTOT > chi2cut. Try to fit without point"
                            << std::endl;
      for (int point = npoints - 1; point >= 0; point--) {
        // Stop if we would remove too many points on the track to meet the minimum hit requirement (or if the track has more
        // than 2 hits removed)
        if ((npoints - removed - 1) < parameters._minClustersOnTrack || removed == 2)
          break;

        // Refit the track without this point
        double newChi2ndofTOT = fitWithoutPoint(*track, point);

        // If the chi2/ndof is significantly better, remove the point permanently CHANGE ME??
        //        if( (chi2ndofTOT - newChi2ndofTOT) > 0 && (chi2ndofTOT - newChi2ndofTOT) > 1. ){
        if ((newChi2ndofTOT - chi2ndofTOT) < chi2ndofTOT) {
          track->remove(point);
          removed++;
          chi2ndofTOT = newChi2ndofTOT;
        }
      }
    }

    streamlog_out(DEBUG8) << "- Track " << candidateTrack.get() << " has " << track->m_clusters.size() << " hits after fit"
                          << std::endl;
    for (unsigned int cluster = 0; cluster < track->m_clusters.size() && streamlog_level(DEBUG8); cluster++) {
      streamlog_out(DEBUG8) << "-- Hit " << cluster << ": [x,y] = [" << track->m_clusters.at(cluster)->getX() << ", "
                            << track->m_clusters.at(cluster)->getY() << "]" << std::endl;
    }

    // Store the track information
    trackContainer.push_back(std::move(track));

    candidateTrack.reset();

  }  // end for candidateTracks

  // Now have all sets of conformal tracks and their chi2/ndof. Decide which tracks to send back, ie. the one with
  // lowest chi2/ndof, and possibly others if they are not clones and have similar chi2 value
  getLowestChi2(finalTracks, trackContainer);
  streamlog_out(DEBUG8) << "getFittedTracks *****" << std::endl;

  // Send back the final set of tracks
  return;
}

// Pick the lowest chi2/ndof KDTrack from a list of possible tracks, and additionally return other tracks in the collection with similar chi2/ndof values that don't share many hits
void ConformalTracking::getLowestChi2(UniqueKDTracks& finalTracks, UniqueKDTracks& trackContainer) {
  streamlog_out(DEBUG8) << "***** getLowestChi2" << std::endl;

  // Get the lowest chi2/ndof value from the given tracks
  //  double lowestChi2ndof = *std::min_element(trackChi2ndofs.begin(),trackChi2ndofs.end());
  UKDTrack& lowestChi2ndofTrack = trackContainer[0];
  double    lowestChi2ndof      = sqrt(lowestChi2ndofTrack->chi2ndof() * lowestChi2ndofTrack->chi2ndof() +
                                       lowestChi2ndofTrack->chi2ndofZS() * lowestChi2ndofTrack->chi2ndofZS());

  for (auto& itTrack : trackContainer) {
    double chi2ndof = sqrt(itTrack->chi2ndof() * itTrack->chi2ndof() + itTrack->chi2ndofZS() * itTrack->chi2ndofZS());

    streamlog_out(DEBUG8) << "- Track " << itTrack.get() << " has " << chi2ndof << " chi2" << std::endl;
    if (chi2ndof < lowestChi2ndof) {
      lowestChi2ndof = itTrack->chi2ndof();
      // lowestChi2ndofTrack = trackContainer[itTrack];
    }
  }

  //finalTracks.push_back(lowestChi2ndofTrack);
  //return;
  streamlog_out(DEBUG8) << "Save more than one track if similar chi2" << std::endl;

  // Loop over all other tracks and decide whether or not to save them
  // Look at the difference in chi2/ndof - we want to keep tracks with similar chi2/ndof. If they
  // are clones then take the longest
  copy_if(std::make_move_iterator(trackContainer.begin()), std::make_move_iterator(trackContainer.end()),
          std::back_inserter(finalTracks), [lowestChi2ndof](UKDTrack const& track) {
            return ((sqrt((track->chi2ndof() * track->chi2ndof() + track->chi2ndofZS() * track->chi2ndofZS())) -
                     lowestChi2ndof) < 10.);
          });

  streamlog_out(DEBUG8) << "Final fitted tracks: " << finalTracks.size() << std::endl;
  if (streamlog_level(DEBUG8)) {
    for (auto& itTrk : finalTracks) {
      streamlog_out(DEBUG8) << "- Track " << itTrk.get() << " has " << itTrk->m_clusters.size() << " hits" << std::endl;
      for (unsigned int cluster = 0; cluster < itTrk->m_clusters.size(); cluster++) {
        streamlog_out(DEBUG8) << "-- Hit " << cluster << ": [x,y] = [" << itTrk->m_clusters.at(cluster)->getX() << ", "
                              << itTrk->m_clusters.at(cluster)->getY() << "]" << std::endl;
      }
    }
  }
  trackContainer.clear();
  streamlog_out(DEBUG8) << "getLowestChi2 *****" << std::endl;

  return;
}

void ConformalTracking::updateCell(SCell const& cell) {
  if (cell->getTo().size() != 0) {
    for (unsigned int i = 0; i < cell->getTo().size(); i++) {
      SCell(cell->getTo()[i])->update(cell);
      updateCell(SCell(cell->getTo()[i]));
    }
  }
}

// Function to extrapolate along a cell in conformal space, producing a fake hit
// a given distance away from the cell endpoint
SKDCluster ConformalTracking::extrapolateCell(SCell const& cell, double distance) {
  // Fake cluster to be returned
  SKDCluster extrapolatedCluster = std::make_shared<KDCluster>();

  // Linear extrapolation of the cell - TODO: check that cell gradients have correct sign and remove checks here
  double gradient = cell->getGradient();
  double deltaU   = sqrt(distance * distance / (1 + gradient * gradient));
  double deltaV   = abs(gradient) * deltaU;

  if ((cell->getStart()->getU() - cell->getEnd()->getU()) > 0)
    deltaU *= (-1.);
  if ((cell->getStart()->getV() - cell->getEnd()->getV()) > 0)
    deltaV *= (-1.);

  extrapolatedCluster->setU(cell->getEnd()->getU() + deltaU);
  extrapolatedCluster->setV(cell->getEnd()->getV() + deltaV);

  return extrapolatedCluster;
}

// NOT CURRENTLY USED - DEPRECATE?
void ConformalTracking::extendTrack(UKDTrack& track, UniqueCellularTracks trackSegments,
                                    std::map<SKDCluster, bool>& /*used*/, std::map<SCell, bool>& /*usedCells*/) {
  //  cout<<"== extending track, have "<<trackSegments.size()<<" candidates"<<endl;
  // For each track segment, perform a kalman filter on the addition points and chose the track extension with the
  // best delta chi2.
  double bestChi2 = 1.e9;
  //double             bestNpoints;
  SharedKDClusters finalHits;
  double           bestChi2Conformal(0.0), bestChi2ZS(0.0);

  //  for(int i=0;i<track->clusters().size();i++)streamlog_out(DEBUG7)<<"- cluster "<<i<<" has u,v "<<track->clusters()[i]->getU()<<","<<track->clusters()[i]->getV()<<std::endl;

  for (auto& trackSegment : trackSegments) {
    //   streamlog_out(DEBUG7)<<"Track extension "<<nTrackExtension<<std::endl;
    SharedKDClusters newHits;

    // Add each of the new clusters
    SKDCluster kdStart = (*trackSegment)[0]->getEnd();
    newHits.push_back(kdStart);
    //   streamlog_out(DEBUG7)<<"- cluster "<<" has u,v "<<kdStart->getU()<<","<<kdStart->getV()<<std::endl;

    for (unsigned int trackCell = 0; trackCell < trackSegment->size() - 2; trackCell++) {
      SKDCluster kdEnd = (*trackSegment)[trackCell]->getStart();
      newHits.push_back(kdEnd);
      //     streamlog_out(DEBUG7)<<"- cluster "<<" has u,v "<<kdEnd->getU()<<","<<kdEnd->getV()<<std::endl;
    }

    double newChi2Conformal, newChi2ZS;
    double newChi2 = fitWithExtension(*track, newHits, newChi2Conformal, newChi2ZS);

    if (newChi2 < bestChi2) {
      finalHits.clear();
      for (auto const& newHit : newHits)
        finalHits.push_back(newHit);
      bestChi2          = newChi2;
      bestChi2Conformal = newChi2Conformal;
      bestChi2ZS        = newChi2ZS;
    }
  }
  if (finalHits.size() == 0)
    return;

  if ((bestChi2ZS - track->chi2ndofZS()) < 100. && (bestChi2Conformal - track->chi2ndof()) < 100.) {
    for (int i = finalHits.size() - 1; i >= 0; i--) {
      //      used[finalHits[i]] = true;
      track->m_clusters.insert(track->m_clusters.begin(), finalHits[i]);
      track->m_nPoints++;
    }
    track->linearRegression();
    track->linearRegressionConformal();
    //   streamlog_out(DEBUG7)<<"- new chi2/ndof: "<<track->chi2ndof()<<", chi2ZS/ndof: "<<track->chi2ndofZS()<<std::endl;
  }
}

//===================================
// Fitting Functions
//===================================

double ConformalTracking::fitWithExtension(KDTrack track, SharedKDClusters hits, double& newChi2, double& newChi2ZS) {
  // Add the point to the track
  for (auto const& hit : hits)
    track.add(hit);

  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal();

  double chi2ndof = track.chi2ndofZS() + track.chi2ndof();
  newChi2         = track.chi2ndof();
  newChi2ZS       = track.chi2ndofZS();
  return chi2ndof;
}

// Add a point to a track and return the delta chi2
void ConformalTracking::fitWithPoint(KDTrack kdTrack, SKDCluster& hit, double& deltaChi2, double& deltaChi2zs) {
  double chi2   = kdTrack.chi2ndof();
  double chi2zs = kdTrack.chi2ndofZS();
  /*if(m_debugPlots){

    double xMeasured = hit->getU();
    double yMeasured = hit->getV();
    double dx = hit->getErrorU();
    double dv = hit->getErrorV();

    if(kdTrack.m_rotated){
      double newxMeasured = yMeasured;
      double newyMeasured = -1. * xMeasured;
      double newdx        = dv;
      double newdv        = dx;
      xMeasured           = newxMeasured;
      yMeasured           = newyMeasured;
      dx                  = newdx;
      dv                  = newdv;
    }
    double residualY = (kdTrack.m_gradient * xMeasured + kdTrack.m_quadratic * xMeasured * xMeasured + kdTrack.m_intercept) - yMeasured;
    // Get the error on the hit position
    double term = kdTrack.m_gradient + 2. * kdTrack.m_quadratic * xMeasured;
    double dy2  = (dv * dv) + (term * term * dx * dx);

    streamlog_out(DEBUG7)<<"- hit has delta chi2 of "<<(residualY * residualY) / (dy2)<<std::endl;

  }*/
  kdTrack.add(hit);
  kdTrack.linearRegression();
  kdTrack.linearRegressionConformal();
  double newchi2   = kdTrack.chi2ndof();
  double newchi2zs = kdTrack.chi2ndofZS();

  deltaChi2   = newchi2 - chi2;
  deltaChi2zs = newchi2zs - chi2zs;
}

// Remove a point from a track and return the delta chi2
double ConformalTracking::fitWithoutPoint(KDTrack track, int point) {
  // Remove the given point from the track
  track.remove(point);

  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal();

  double chi2ndof   = track.chi2ndof();
  double chi2ndofZS = track.chi2ndofZS();

  return chi2ndof + chi2ndofZS;
}

//===================================
// Utility Functions
//===================================

// Get a collection from the event object
void ConformalTracking::getCollection(LCCollection*& collection, std::string const& collectionName, LCEvent* evt) {
  try {
    collection = evt->getCollection(collectionName);
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG7) << "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

// Function to check if two KDtracks contain several hits in common
int ConformalTracking::overlappingHits(const UKDTrack& track1, const UKDTrack& track2) {
  int nHitsInCommon = 0;
  for (size_t hit = 0; hit < track1->m_clusters.size(); hit++) {
    if (std::find(track2->m_clusters.begin(), track2->m_clusters.end(), track1->m_clusters[hit]) != track2->m_clusters.end())
      nHitsInCommon++;
  }
  return nHitsInCommon;
}

//===================================
// Debug Functions for Reconstruction
//===================================

// Check if cellular tracks were valid or whether any contain unallowed cell combinations
void ConformalTracking::checkUnallowedTracks(UniqueCellularTracks candidateTracks, Parameters const& parameters) {
  double nTracks          = candidateTracks.size();
  double nUnallowedTracks = 0.;

  for (unsigned int itTrack = 0; itTrack < nTracks; itTrack++) {
    const auto&  track  = candidateTracks[itTrack];
    unsigned int nCells = track->size();

    for (unsigned int trackCell = 0; trackCell < nCells - 1; trackCell++) {
      double cellAngle   = (*track)[trackCell]->getAngle((*track)[trackCell + 1]);
      double cellAngleRZ = (*track)[trackCell]->getAngleRZ((*track)[trackCell + 1]);

      if (cellAngle > parameters._maxCellAngle || cellAngleRZ > parameters._maxCellAngleRZ) {
        nUnallowedTracks++;
        break;
      }
    }
  }

  streamlog_out(DEBUG7) << " Produced " << nTracks << " tracks, " << nUnallowedTracks
                        << " of which were bad combinations: " << nUnallowedTracks / nTracks << " %" << std::endl;
}

// Debug function - checks if a track will be associated to an MC particle or not
double ConformalTracking::checkReal(UKDTrack& track, std::map<MCParticle*, bool>& reconstructed,
                                    std::map<MCParticle*, SharedKDClusters> MCparticleHits) {
  // Store all mcparticles associated to this track
  std::vector<MCParticle*>      particles;
  std::map<MCParticle*, double> particleHits;
  double                        nHits = 0.;

  // Get the clusters from this track
  SharedKDClusters clusters = track->m_clusters;

  // Loop over all hits and see which particle they are associated to
  for (auto const& cluster : clusters) {
    // Get the hit
    nHits++;

    // If we already have hits on this particle, then just increment the counter
    if (particleHits.count(kdParticles[cluster]))
      particleHits[kdParticles[cluster]]++;
    else {
      // Otherwise register the new particle and set the number of hits to 1
      particles.push_back(kdParticles[cluster]);
      particleHits[kdParticles[cluster]] = 1;
    }
  }

  // Now look how many hits are on each particle and calculate the purity
  double      bestHits     = 0.;
  MCParticle* bestParticle = nullptr;
  for (auto& particle : particles) {
    if (particleHits[particle] > bestHits) {
      bestHits     = particleHits[particle];
      bestParticle = particle;
    }
  }

  // Calculate the purity
  double purity = bestHits / nHits;
  streamlog_out(DEBUG9) << "Number of hits on track: " << nHits << ". Good hits: " << bestHits << ". Purity: " << purity
                        << ". Pt: "
                        << sqrt(bestParticle->getMomentum()[0] * bestParticle->getMomentum()[0] +
                                bestParticle->getMomentum()[1] * bestParticle->getMomentum()[1])
                        << ". Pt estimate: " << track->pt() << ". Track chi2/ndof: " << track->chi2ndof()
                        << ". Chi2/ndof in SZ fit: " << track->chi2ndofZS() << std::endl;

  // Check the track pt estimate
  TLorentzVector mc_helper;
  mc_helper.SetPxPyPzE(bestParticle->getMomentum()[0], bestParticle->getMomentum()[1], bestParticle->getMomentum()[2],
                       bestParticle->getEnergy());

  streamlog_out(DEBUG7) << "Particle at theta = " << mc_helper.Theta() * 180. / M_PI << ", with pt = " << mc_helper.Pt()
                        << ". pt estimate: " << track->pt() << std::endl;

  // Check if any hits are missing
  SharedKDClusters mcHits     = MCparticleHits[bestParticle];
  int              uniqueHits = getUniqueHits(mcHits);

  streamlog_out(DEBUG7) << "Track should contain " << uniqueHits << " hits, "
                        << ((uniqueHits > bestHits) ? "is missing hits" : "all hits found.") << std::endl;

  SharedKDClusters trackHits = track->m_clusters;
  for (size_t th = 0; th < trackHits.size(); th++) {
    streamlog_out(DEBUG9) << "Hit " << th << " u = " << trackHits[th]->getU() << " v = " << trackHits[th]->getV()
                          << " x = " << trackHits[th]->getX() << " y = " << trackHits[th]->getY()
                          << " z = " << trackHits[th]->getZ() << std::endl;

    // Check which particle the hit belongs to
    MCParticle* particle  = kdParticles[trackHits[th]];
    double      mcVertexX = particle->getVertex()[0];
    double      mcVertexY = particle->getVertex()[1];
    double      mcVertexR = sqrt(pow(mcVertexX, 2) + pow(mcVertexY, 2));

    streamlog_out(DEBUG9) << "Come from particle with id " << particle->getPDG() << " produced at vertexR = " << mcVertexR;

    if (particle->getParents().size() != 0) {
      streamlog_out(DEBUG9) << ". Comes from particle with id " << particle->getParents()[0]->getPDG();
    }
    streamlog_out(DEBUG7) << std::endl;
  }

  //  for(int itCluster=0;itCluster<clusters.size();itCluster++)streamlog_out(DEBUG7)<<"Hit "<<itCluster<<" has position: "<<clusters[itCluster]->getU()<<","<<clusters[itCluster]->getV()<<std::endl;
  streamlog_out(DEBUG7) << "== Terms in conformal fit: " << track->m_intercept << ", " << track->m_gradient << ", "
                        << track->m_quadratic << std::endl;
  streamlog_out(DEBUG7) << "== Terms in zs fit: " << track->m_interceptZS << ", " << track->m_gradientZS << std::endl;

  track->linearRegressionConformal(true);
  track->calculateChi2SZ(nullptr, true);
  if (purity >= m_purity) {
    reconstructed[bestParticle] = true;
  }

  // Return the purity
  return purity;
}

// Debug function - gets number of unique hits
int ConformalTracking::getUniqueHits(SharedKDClusters hits) {
  int nUniqueHits = 0;

  for (size_t iHit = 0; iHit < hits.size(); iHit++) {
    // Get the cluster
    SKDCluster hit    = hits[iHit];
    bool       sameID = false;

    // Check if any following clusters have the same detector id
    for (size_t iHitLater = iHit + 1; iHitLater < hits.size(); iHitLater++) {
      // Get the new cluster
      SKDCluster nhit = hits[iHitLater];

      if (hit->sameLayer(nhit))
        sameID = true;
    }

    // Didn't find another hit with the same ID
    if (!sameID)
      nUniqueHits++;
  }

  return nUniqueHits;
}

// Make the track in the same way as the pattern recognition would, but
// without the inclusion of unassociated hits. See which criteria fail
void ConformalTracking::checkReconstructionFailure(MCParticle*                             particle,
                                                   std::map<MCParticle*, SharedKDClusters> particleHits,
                                                   UKDTree& nearestNeighbours, Parameters const& parameters) {
  // Get the hits for this MC particle
  SharedKDClusters trackHits = particleHits[particle];

  // Sort the hits from larger to smaller radius
  std::sort(trackHits.begin(), trackHits.end(), sort_by_radiusKD);
  streamlog_out(DEBUG7) << "Track contains " << getUniqueHits(trackHits) << " unique hits. Track size is "
                        << trackHits.size() << std::endl;

  for (size_t th = 0; th < trackHits.size(); th++)
    streamlog_out(DEBUG7) << "Hit " << th << " u = " << trackHits[th]->getU() << " v = " << trackHits[th]->getV()
                          << " x = " << trackHits[th]->getX() << " y = " << trackHits[th]->getY()
                          << " z = " << trackHits[th]->getZ() << ". " << (trackHits[th]->used() ? "Used!" : "") << std::endl;

  // Take the seed hit and build cell to 2nd
  SKDCluster       seedHit = trackHits[0];
  SharedKDClusters results;
  streamlog_out(DEBUG7) << "- getting nearest neighbours for seed" << std::endl;
  //  nearestNeighbours->allNeighboursInRadius(seedHit, parameters._maxDistance, results);
  double theta = seedHit->getTheta();
  nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);
  streamlog_out(DEBUG7) << "- got nearest neighbours for seed" << std::endl;

  // If we don't produce the seed cell we have failed
  if (std::find(results.begin(), results.end(), trackHits[1]) == results.end()) {
    streamlog_out(DEBUG7) << "- seed cell not produced, neighbour not inside search window" << std::endl;
    double deltaU = fabs(trackHits[1]->getU() - trackHits[0]->getU());
    double deltaV = fabs(trackHits[1]->getV() - trackHits[0]->getV());
    streamlog_out(DEBUG7) << "- distance between hits is " << sqrt(deltaU * deltaU + deltaV * deltaV)
                          << " and search distance is " << parameters._maxDistance << std::endl;
    streamlog_out(DEBUG7) << "- theta of hits is " << trackHits[0]->getTheta() << " and " << trackHits[1]->getTheta()
                          << ", delta theta = " << trackHits[0]->getTheta() - trackHits[1]->getTheta() << std::endl;
    //    return;
  }

  if (trackHits[1]->getR() >= trackHits[0]->getR()) {
    streamlog_out(DEBUG7) << "- seed cell not produced, neighbour at higher radius" << std::endl;
    //    return;
  }

  double length = sqrt((trackHits[1]->getU() - trackHits[0]->getU()) * (trackHits[1]->getU() - trackHits[0]->getU()) +
                       (trackHits[1]->getV() - trackHits[0]->getV()) * (trackHits[1]->getV() - trackHits[0]->getV()));
  if (length > parameters._maxDistance) {
    streamlog_out(DEBUG7) << "- seed cell not produced, neighbour too far away" << std::endl;
    //    return;
  }

  // Make the seed cell and extrapolate until all hits found
  auto        seedCell = std::make_shared<Cell>(trackHits[0], trackHits[1]);
  SharedCells cells;
  cells.push_back(seedCell);
  streamlog_out(DEBUG7) << "have seed cell" << std::endl;
  // CHECK CELL LENGTH!! TO DO

  int hitNumber = 1;  // Current hit after the seed cell has been formed
  //bool trackBuilding = true;
  while (hitNumber < (int(trackHits.size()) - 1)) {
    // Extend the current cell
    double searchDistance = parameters._maxDistance;  //hit->getR();
    //    if(searchDistance > trackHits[hitNumber]->getR()) searchDistance = trackHits[hitNumber]->getR();

    // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
    SKDCluster const& fakeHit =
        extrapolateCell(cells.back(), searchDistance / 2.);  // TODO: make this search a function of radius
    SharedKDClusters results2;
    nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25 * searchDistance / 2., results2);

    // Check if we found the hit we wanted
    if (std::find(results2.begin(), results2.end(), trackHits[hitNumber + 1]) == results2.end()) {
      streamlog_out(DEBUG7) << "- could not find hit number " << hitNumber + 1 << " inside the search window" << std::endl;
      double deltaU = fabs(trackHits[hitNumber + 1]->getU() - trackHits[hitNumber]->getU());
      double deltaV = fabs(trackHits[hitNumber + 1]->getV() - trackHits[hitNumber]->getV());
      streamlog_out(DEBUG7) << "- distance between hits is " << sqrt(deltaU * deltaU + deltaV * deltaV)
                            << " and search distance is " << parameters._maxDistance << std::endl;
      //      for(int c=0;c<cells.size();c++) delete cells[c];
      //      return;
    }

    // Check radial conditions
    if (trackHits[hitNumber + 1]->getR() >= trackHits[hitNumber]->getR()) {
      streamlog_out(DEBUG7) << "- cell " << hitNumber << " not produced, neighbour at higher radius" << std::endl;
      //      return;
    }

    // Make the cell for extrapolation in the next round
    auto cell = std::make_shared<Cell>(trackHits[hitNumber], trackHits[hitNumber + 1]);
    streamlog_out(DEBUG7) << "have new cell" << std::endl;

    // Check if the cell would be killed by cuts
    if (cells.back()->getAngle(cell) > parameters._maxCellAngle ||
        cells.back()->getAngleRZ(cell) > parameters._maxCellAngleRZ) {
      if (cells.back()->getAngle(cell) > parameters._maxCellAngle)
        streamlog_out(DEBUG7) << "- cell " << hitNumber << " killed by angular cut. Angle is "
                              << cells.back()->getAngle(cell) << std::endl;
      if (cells.back()->getAngleRZ(cell) > parameters._maxCellAngleRZ)
        streamlog_out(DEBUG7) << "- cell " << hitNumber << " killed by RZ angular cut. Angle is "
                              << cells.back()->getAngleRZ(cell) << std::endl;
      //      for(int c=0;c<cells.size();c++) delete cells[c];
      //      return;
    }

    // Set the information about which cell this new cell is attached to and store it
    cell->setFrom(SCell(cells.back()));
    cells.back()->setTo(cell);
    //    cell->setWeight(hitNumber);
    cells.push_back(cell);

    // Check if we have finished
    //    if(hitNumber == (trackHits.size()-1)) trackBuilding = false;
    hitNumber++;
  }

  // Now test the built-in functions to make CellularTracks from cells
  UniqueCellularTracks  cellularTracks;
  std::map<SCell, bool> usedCells;
  createTracksNew(cellularTracks, cells.back(), usedCells);

  if (cellularTracks.size() != 1) {
    streamlog_out(DEBUG7) << "- cellular track not produced from cells. Returned " << cellularTracks.size() << " candidates"
                          << std::endl;
    return;
  }

  // Check that all of the hits are on the track
  if (cellularTracks[0]->size() != (trackHits.size() - 1)) {
    streamlog_out(DEBUG7) << "- failed to put all cells on track. Cellular track size: " << cellularTracks[0]->size()
                          << std::endl;
    streamlog_out(DEBUG7) << "- number of cells held: " << cells.size() << std::endl;
    return;
  }

  // Now test the built-in functions to make KDTracks from CellularTracks
  UniqueKDTracks finalTracks;
  getFittedTracks(finalTracks, cellularTracks, usedCells, parameters);
  cellularTracks.clear();  // no longer used

  if (finalTracks.size() != 1) {
    streamlog_out(DEBUG7) << "- kd track not produced from cellular track. Returned " << finalTracks.size() << " candidates"
                          << std::endl;
    return;
  }

  // Test the chi2 criteria
  UKDTrack& mcTrack = finalTracks[0];
  mcTrack->linearRegression(false);
  mcTrack->linearRegressionConformal(true);

  streamlog_out(DEBUG7) << "== Track chi2/ndof is " << mcTrack->chi2ndof() << ", ZS chi2/ndof is " << mcTrack->chi2ndofZS()
                        << std::endl;
  streamlog_out(DEBUG7) << "== Terms in conformal fit: " << mcTrack->m_intercept << ", " << mcTrack->m_gradient << ", "
                        << mcTrack->m_quadratic << std::endl;
  streamlog_out(DEBUG7) << "== Terms in zs fit: " << mcTrack->m_interceptZS << ", " << mcTrack->m_gradientZS << std::endl;
  mcTrack->calculateChi2SZ(nullptr, true);

  if (mcTrack->chi2ndof() > parameters._chi2cut || mcTrack->chi2ndofZS() > parameters._chi2cut) {
    if (mcTrack->chi2ndof() > parameters._chi2cut)
      streamlog_out(DEBUG7) << "- track killed by chi2/ndof cut. Track chi2/ndof is " << mcTrack->chi2ndof() << std::endl;
    if (mcTrack->chi2ndofZS() > parameters._chi2cut)
      streamlog_out(DEBUG7) << "- track killed by ZS chi2/ndof cut. Track chi2/ndof in ZS is " << mcTrack->chi2ndofZS()
                            << std::endl;
    return;
  }

  streamlog_out(DEBUG7) << "== Should have produced this track" << std::endl;
  return;
}

// Draw a line on the current canvas
void ConformalTracking::drawline(SKDCluster const& hitStart, SKDCluster const& hitEnd, int colour, int style) {
  TLine* line = new TLine(hitStart->getU(), hitStart->getV(), hitEnd->getU(), hitEnd->getV());
  line->SetLineColor(colour);
  line->SetLineStyle(style);
  line->Draw();
}

void ConformalTracking::fillCollectionIndexVectors() {
  std::vector<std::string> vertexCombinedHitCollections = m_inputVertexBarrelCollections;
  vertexCombinedHitCollections.insert(vertexCombinedHitCollections.end(), m_inputVertexEndcapCollections.begin(),
                                      m_inputVertexEndcapCollections.end());

  std::vector<std::pair<std::vector<int>*, std::vector<std::string>*>> indexCollectionMap = {
      std::make_pair(&m_allHits, &m_inputTrackerHitCollections),
      std::make_pair(&m_trackerHits, &m_inputMainTrackerHitCollections),
      std::make_pair(&m_vertexBarrelHits, &m_inputVertexBarrelCollections),
      std::make_pair(&m_vertexEndcapHits, &m_inputVertexEndcapCollections),
      std::make_pair(&m_vertexCombinedHits, &vertexCombinedHitCollections)};

  for (auto& pairIndexCol : indexCollectionMap) {
    auto* indices = pairIndexCol.first;
    for (auto const& colName : *(pairIndexCol.second)) {
      auto it = std::find(m_inputTrackerHitCollections.begin(), m_inputTrackerHitCollections.end(), colName);
      if (it == m_inputTrackerHitCollections.end()) {
        std::stringstream error;
        error << this->name() << ": Collection name \"" << colName
              << "\" not found in TrackerHitCollectionNames parameter! Check your config";
        throw marlin::ParseException(error.str());
      }
      indices->push_back(it - m_inputTrackerHitCollections.begin());
    }
  }

  for (auto& pairIndexCol : indexCollectionMap) {
    auto*             indices  = pairIndexCol.first;
    auto*             colNames = pairIndexCol.second;
    std::stringstream indexList;
    for (size_t i = 0; i < indices->size(); ++i) {
      indexList << colNames->at(i) << "=" << indices->at(i) << "   ";
    }
    streamlog_out(MESSAGE) << "IndexList " << indexList.str() << std::endl;
    ;
  }
}

void ConformalTracking::runStep(SharedKDClusters& kdClusters, UKDTree& nearestNeighbours, UniqueKDTracks& conformalTracks,
                                std::map<int, SharedKDClusters> const& collectionClusters, Parameters const& parameters) {
  auto stopwatch = TStopwatch();
  stopwatch.Start(false);

  if (parameters._combine) {
    combineCollections(kdClusters, nearestNeighbours, parameters._collections, collectionClusters);
  }

  if (parameters._build) {
    bool       caughtException = false;
    Parameters thisParameters(parameters);
    do {
      caughtException = false;
      try {
        buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, thisParameters, parameters._radialSearch,
                       parameters._vertexToTracker);
      } catch (TooManyTracksException& e) {
        if (not m_retryTooManyTracks || thisParameters._tightenStep >= 10) {
          streamlog_out(ERROR) << "Skipping event" << std::endl;
          throw;
        }
        streamlog_out(MESSAGE) << "caught too many tracks, tightening parameters" << std::endl;
        thisParameters.tighten();
        caughtException = true;
      }
    } while (caughtException);

    stopwatch.Stop();
    if (m_debugTime)
      streamlog_out(DEBUG9) << "Step " << parameters._step << " buildNewTracks took " << stopwatch.RealTime() << " seconds"
                            << std::endl;
    stopwatch.Reset();
  }
  if (parameters._extend) {
    extendTracks(conformalTracks, kdClusters, nearestNeighbours, parameters);
    stopwatch.Stop();
    if (m_debugTime)
      streamlog_out(DEBUG9) << "Step " << parameters._step << " extendTracks took " << stopwatch.RealTime() << " seconds"
                            << std::endl;
    stopwatch.Reset();
  }

  if (parameters._sortTracks) {
    // FIXME: Still needed?
    // Sort by pt (low to hight)
    std::sort(conformalTracks.begin(), conformalTracks.end(), sort_by_pt);
  }

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }
}
