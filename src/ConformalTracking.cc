#include "ConformalTracking.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/BitSet32.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>
#include "UTIL/LCTrackerConf.h"

#include <marlinutil/GeometryUtil.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <marlin/Exceptions.h>
#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"

#include "CLHEP/Vector/TwoVector.h"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLinearFitter.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include "TStopwatch.h"

#include "Cell.h"
#include "KDTree.h"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace dd4hep;
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

ConformalTracking::ConformalTracking() : Processor("ConformalTracking") {
  // Processor description
  _description = "ConformalTracking constructs tracks using a combined conformal mapping and cellular automaton approach.";

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
  registerProcessorParameter("ThetaRange", "Angular range for initial cell seeding", m_thetaRange, double(0.1));
  registerProcessorParameter("MaxCellAngle", "Cut on angle between two cells for cell to be valid", m_maxCellAngle,
                             double(0.035));
  registerProcessorParameter("MaxCellAngleRZ", "Cut on angle between two cells in RZ for cell to be valid", m_maxCellAngleRZ,
                             double(0.035));
  registerProcessorParameter("MaxDistance", "Maximum length of a cell (max. distance between two hits)", m_maxDistance,
                             double(0.015));
  registerProcessorParameter("MaxChi2", "Maximum chi2/ndof for linear conformal tracks", m_chi2cut, double(300.));
  registerProcessorParameter("MinClustersOnTrack", "Minimum number of clusters to create a track", m_minClustersOnTrack,
                             int(6));
  registerProcessorParameter("trackPurity", "Purity value used for checking if tracks are real or not", m_purity,
                             double(0.75));
  registerProcessorParameter("MaxChi2Increase", "Chi2 increase when adding new hits to a track", m_chi2increase,
                             double(10.));
}

void ConformalTracking::init() {
  // Print the initial parameters
  printParameters();

  // Reset counters
  m_runNumber   = 0;
  m_eventNumber = 0;

  // Set up the track fit factory
  const gear::GearMgr* fakeGear = 0;
  trackFactory                  = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", fakeGear, "");
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
  debugSeed = NULL;

  // Initialise histograms (if debug plotting on)
  if (m_debugPlots) {
    // Automatically add histograms to output root file
    AIDAProcessor::histogramFactory(this);

    // Histogram initailisation
    m_szDistribution  = new TH2F("m_szDistribution", "m_szDistribution", 20000, -100, 100, 200, -10, 10);
    m_uvDistribution  = new TH2F("m_uvDistribution", "m_uvDistribution", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    m_xyDistribution  = new TH2F("m_xyDistribution", "m_xyDistribution", 500, -1500, 1500, 500, -1500, 1500);
    m_xyzDistribution = new TH3F("m_xyzDistribution", "m_xyzDistribution", 50, 0, 100, 50, 0, 100, 100, 0, 25);

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

  streamlog_out(DEBUG7) << "Event number: " << m_eventNumber << std::endl;

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
    streamlog_out(DEBUG7) << "Collection " << m_inputTrackerHitCollections[collection] << " contains "
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
  std::map<SKDCluster, MCParticle*>       kdParticles;    // Link from conformal hit to MC particle
  std::map<MCParticle*, SharedKDClusters> particleHits;   // List of conformal hits on each MC particle
  std::map<MCParticle*, bool>             reconstructed;  // Check for MC particles
  SharedKDClusters debugHits;                             // Debug hits for plotting
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
      kdhit->setDetectorInfo(subdet, side, layer);

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
      auto mcTrack = std::unique_ptr<KDTrack>(new KDTrack());
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
  auto             stopwatch         = std::unique_ptr<TStopwatch>(new TStopwatch());

  // Build tracks in the vertex barrel

  //Create ParameterStruct
  Parameters initialParameters(m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut, m_minClustersOnTrack, m_maxDistance, true,
                               false);

  stopwatch->Start(false);
  std::vector<int> vertexHits = {0};
  combineCollections(kdClusters, nearestNeighbours, vertexHits, collectionClusters);
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, initialParameters);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters) {
      thisCluster->used(true);
    }
  }

  stopwatch->Stop();
  streamlog_out(DEBUG7) << "Building vertex barrel tracks took " << stopwatch->RealTime() << " seconds" << std::endl;
  stopwatch->Reset();

  // Extend through the endcap
  stopwatch->Start(false);
  std::vector<int> vertexEndcapHits = {1};
  combineCollections(kdClusters, nearestNeighbours, vertexEndcapHits, collectionClusters);
  extendTracks(conformalTracks, kdClusters, nearestNeighbours, initialParameters);
  stopwatch->Stop();

  streamlog_out(DEBUG7) << "Extending through vertex endcap took " << stopwatch->RealTime() << " seconds" << std::endl;
  stopwatch->Reset();

  // Make combined vertex tracks
  stopwatch->Start(false);
  std::vector<int> vertexCombinedHits = {0, 1};
  combineCollections(kdClusters, nearestNeighbours, vertexCombinedHits, collectionClusters);
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, initialParameters);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }
  stopwatch->Stop();
  streamlog_out(DEBUG7) << "Building vertex tracks took " << stopwatch->RealTime() << " seconds" << std::endl;
  stopwatch->Reset();

  // Make leftover tracks in the vertex with lower requirements

  // Lower cell angle criteria
  //m_highPTfit = false;

  Parameters lowerCellAngleParameters(m_maxCellAngle * 5.0, m_maxCellAngleRZ * 5.0, m_chi2cut, m_minClustersOnTrack,
                                      m_maxDistance, true, false);
  //  m_chi2cut*=20.;
  stopwatch->Start(false);
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, lowerCellAngleParameters, true);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }

  stopwatch->Stop();

  streamlog_out(DEBUG7) << "Building low pt vertex tracks (1) took " << stopwatch->RealTime() << " seconds" << std::endl;
  stopwatch->Reset();
  stopwatch->Start(false);

  lowerCellAngleParameters._maxCellAngle *= 2.;
  lowerCellAngleParameters._maxCellAngleRZ *= 2.;

  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, lowerCellAngleParameters, true);
  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }

  stopwatch->Stop();

  streamlog_out(DEBUG7) << "Building low pt vertex tracks (2) took " << stopwatch->RealTime() << " seconds" << std::endl;
  stopwatch->Reset();

  // Lower number of hits on track

  Parameters lowNumberHitsParameters(m_maxCellAngle * 10, m_maxCellAngleRZ * 10, m_chi2cut * 20.,
                                     /*minClustersOnTrack=*/4, m_maxDistance, true, false);

  stopwatch->Start(false);
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, lowNumberHitsParameters, true);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }
  stopwatch->Stop();

  streamlog_out(DEBUG7) << "Building low pt vertex tracks with 4 hits took " << stopwatch->RealTime() << " seconds"
                        << std::endl;
  stopwatch->Reset();

  // Sort by pt (low to hight)
  std::sort(conformalTracks.begin(), conformalTracks.end(), sort_by_pt);

  //m_highPTfit = false;

  lowNumberHitsParameters._chi2cut = m_chi2cut;

  // Extend them through the inner and outer trackers
  std::vector<int> trackerHits = {2, 3, 4, 5};
  stopwatch->Start(false);
  combineCollections(kdClusters, nearestNeighbours, trackerHits, collectionClusters);
  extendTracks(conformalTracks, kdClusters, nearestNeighbours, lowNumberHitsParameters);
  stopwatch->Stop();

  streamlog_out(DEBUG7) << "Extending through trackers took " << stopwatch->RealTime() << " seconds" << std::endl;
  stopwatch->Reset();

  //  extendHighPT(conformalTracks, kdClusters, nearestNeighbours);
  // Mark hits from "good" tracks as being used
  //  for (unsigned int itTrack = 0; itTrack < conformalTracks.size(); itTrack++) {
  //    for (unsigned int itHit = 0; itHit < conformalTracks[itTrack]->m_clusters.size(); itHit++)
  //      conformalTracks[itTrack]->m_clusters[itHit]->used(true);
  //  }
  //  combineCollections(kdClusters, nearestNeighbours, trackerHits, collectionClusters);

  // Finally reconstruct displaced tracks

  std::vector<int> allHits = {0, 1, 2, 3, 4, 5};
  stopwatch->Start(false);
  combineCollections(kdClusters, nearestNeighbours, allHits, collectionClusters);
  double factor          = 10;
  bool   caughtException = false;
  do {
    caughtException = false;
    try {
      Parameters displacedParameters(m_maxCellAngle * factor, m_maxCellAngleRZ * factor, m_chi2cut * factor, 5, 0.015, false,
                                     true);
      streamlog_out(DEBUG8) << "Building new tracks displaced: " << factor << std::endl;
      buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, displacedParameters, true, false);
    } catch (TooManyTracksException& e) {
      streamlog_out(MESSAGE) << "caught too many tracks, tightening parameters" << std::endl;
      caughtException = true;
      factor -= 1;
      if (not m_retryTooManyTracks && factor <= 0) {
        streamlog_out(ERROR) << "Skipping event" << std::endl;
        throw;
      }
    }
  } while (caughtException);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }
  stopwatch->Stop();
  streamlog_out(DEBUG7) << "Building displaced tracks with all detectors took " << stopwatch->RealTime() << " seconds"
                        << std::endl;
  stopwatch->Reset();

  //*/

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
  streamlog_out(DEBUG7) << "*** CA has made " << conformalTracksFinal.size()
                        << (conformalTracksFinal.size() == 1 ? " track ***" : " tracks ***") << std::endl;

  // Loop over all track candidates
  for (auto& conformalTrack : conformalTracksFinal) {
    // Vector of all the hits on the track
    streamlog_out(DEBUG5) << "Made a track with " << conformalTrack->m_clusters.size() << " hits" << std::endl;

    // Make the LCIO track hit vector
    EVENT::TrackerHitVec trackHits;
    for (auto const& cluster : conformalTrack->m_clusters) {
      trackHits.push_back(kdClusterMap[cluster]);
    }
    // Add kalman filtered hits
    if (conformalTrack->kalmanTrack() != NULL) {
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

    // Try to fit
    int fitError =
        MarlinTrk::createFinalisedLCIOTrack(marlinTrack.get(), trackHits, track.get(), MarlinTrk::IMarlinTrack::forward,
                                            covMatrix, m_magneticField, m_maxChi2perHit);

    // Check track quality - if fit fails chi2 will be 0. For the moment add hits by hand to any track that fails the track fit, and store it as if it were ok...
    if (track->getChi2() == 0.) {
      streamlog_out(DEBUG7) << "Fit failed. Track has " << track->getTrackerHits().size() << " hits" << std::endl;
      streamlog_out(DEBUG7) << "Fit fail error " << fitError << std::endl;
      continue;
    }  //*/

    /*for (unsigned int p = 0; p < trackHits.size(); p++) {
      track->addHit(trackHits[p]);
    }  //*/

    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2] = trackHits.size();

    // calculate purities and check if track has been reconstructed
    if (m_debugPlots) {
      m_conformalChi2->Fill(conformalTrack->chi2ndof());
      streamlog_out(DEBUG7) << "-------------------- New TRACK --------------------" << std::endl;
      //streamlog_out(DEBUG7) << " LCIO track fit chi2 is "<<track->getChi2()<<std::endl;
      double purity = checkReal(conformalTrack, kdParticles, reconstructed, particleHits);
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
      streamlog_out(DEBUG7) << "-------------------- New PARTICLE --------------------" << std::endl;
      // List the pt
      double particlePt = sqrt(mcParticle->getMomentum()[0] * mcParticle->getMomentum()[0] +
                               mcParticle->getMomentum()[1] * mcParticle->getMomentum()[1]);
      streamlog_out(DEBUG7) << "Particle pt: " << particlePt << std::endl;
      checkReconstructionFailure(mcParticle, particleHits, nearestNeighbours_debug, initialParameters);
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
    streamlog_out(DEBUG7) << "Reconstructed " << nReconstructed << " particles out of " << nReconstructed + nUnreconstructed
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

  //FIXME trackFactory is leaking Memory, but probably a MarlinTRK issue
}

//===================================
// Tracking strategies
//===================================

// Combine collections
void ConformalTracking::combineCollections(SharedKDClusters& kdClusters, UKDTree& nearestNeighbours,
                                           std::vector<int> const& combination,
                                           std::map<int, SharedKDClusters> const& collectionClusters) {
  // Clear the input objects
  kdClusters.clear();
  nearestNeighbours.reset(nullptr);

  // Loop over all given collections
  for (unsigned int i = 0; i < combination.size(); i++) {
    // Copy the clusters to the output vector
    const SharedKDClusters& clusters = collectionClusters.at(combination[i]);
    int                     nhits    = clusters.size();
    for (int hit = 0; hit < nhits; hit++) {
      kdClusters.push_back(clusters[hit]);
    }
  }

  streamlog_out(DEBUG7) << "Combined collection has " << kdClusters.size() << " hits" << std::endl;

  // Sort the KDClusters from larger to smaller radius
  std::sort(kdClusters.begin(), kdClusters.end(), sort_by_radiusKD);

  // Make the binary search tree. This tree class contains two binary trees - one sorted by u-v and the other by theta
  nearestNeighbours = UKDTree(new KDTree(kdClusters, m_thetaRange, m_sortTreeResults));
}

// Take a collection of hits and try to produce tracks out of them
void ConformalTracking::buildNewTracks(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                       UKDTree& nearestNeighbours, Parameters const& parameters, bool radialSearch,
                                       bool vertexToTracker) {
  streamlog_out(DEBUG7) << "BUILDING new tracks" << std::endl;

  // Sort the input collection by radius - higher to lower if starting with the vertex detector (high R in conformal space)
  std::sort(collection.begin(), collection.end(), (vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD));

  // Loop over all hits, using each as a seed to produce a new track
  unsigned int nKDHits = collection.size();
  for (unsigned int nKDHit = 0; nKDHit < nKDHits; nKDHit++) {
    // Get the kdHit and check if it has already been used (assigned to a track)
    SKDCluster kdhit = collection[nKDHit];
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "Starting to seed with debug cluster" << std::endl;
    if (kdhit->used())
      continue;
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
    if (radialSearch)
      nearestNeighbours->allNeighboursInRadius(kdhit, parameters._maxDistance, results);
    else
      nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);

    // Sort the neighbours by radius
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "- picked up " << results.size() << " neighbours from " << (radialSearch ? "radial" : "theta")
                            << " search" << std::endl;
    if (results.size() == 0)
      continue;
    std::sort(results.begin(), results.end(), (vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD));

    // Objects to hold cells
    SharedCells cells;

    // Make seed cells pointing inwards/outwards (conformal space)
    for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
      // Get the neighbouring hit
      SKDCluster const& nhit = results[neighbour];

      // Check that it is not used, is not on the same detector layer, is not in the opposite side of the detector and points inwards
      if (nhit->used()) {
        if (debugSeed && kdhit == debugSeed)
          streamlog_out(DEBUG7) << "- used" << std::endl;
        continue;
      }
      if (kdhit->sameLayer(nhit)) {
        if (debugSeed && kdhit == debugSeed)
          streamlog_out(DEBUG7) << "- same layer" << std::endl;
        continue;
      }
      if (nhit->endcap() && kdhit->endcap() && (nhit->forward() != kdhit->forward())) {
        if (debugSeed && kdhit == debugSeed)
          streamlog_out(DEBUG7) << "- not pointing in the same direction" << std::endl;
        continue;
      }
      if ((vertexToTracker && nhit->getR() >= kdhit->getR()) || (!vertexToTracker && nhit->getR() <= kdhit->getR())) {
        if (debugSeed && kdhit == debugSeed)
          streamlog_out(DEBUG7) << "- radial conditions not met" << std::endl;
        continue;
      }

      // Check if the cell would be too long (hit very far away)
      double length2 = ((kdhit->getU() - nhit->getU()) * (kdhit->getU() - nhit->getU()) +
                        (kdhit->getV() - nhit->getV()) * (kdhit->getV() - nhit->getV()));
      if (length2 > parameters._maxDistance * parameters._maxDistance) {
        continue;
      }
      // Create the new seed cell
      cells.emplace_back(std::make_shared<Cell>(kdhit, nhit));
      auto const& cell = cells.back();

      //      if (cell->doca() > 0.01) {
      //        continue;
      //      }

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

    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG7) << "- after extension, have " << cells.size() << " cells" << std::endl;

    // Now have all cells stemming from this seed hit. If it is possible to produce a track (ie. cells with depth X) then we will now...
    //      if(depth < (m_minClustersOnTrack-1)) continue; // TODO: check if this is correct

    // We create all acceptable tracks by looping over all cells with high enough weight to create
    // a track and trace their route back to the seed hit. We then have to choose the best candidate
    // at the end (by minimum chi2 of a linear fit)
    std::map<SCell, bool> usedCells;
    std::map<SCell, bool> usedCells2;
    UniqueKDTracks cellTracks;

    // Sort Cells from highest to lowest weight
    std::sort(cells.begin(), cells.end(), sort_by_cellWeight);

    // Create tracks by following a path along cells
    int nCells = cells.size();
    for (int itCell = 0; itCell < nCells; itCell++) {
      // Check if this cell has already been used
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG7) << "-- looking at cell " << itCell << std::endl;
      if (usedCells.count(cells[itCell]))
        continue;

      // Check if this cell could produce a track (is on a long enough chain)
      if (cells[itCell]->getWeight() < (parameters._minClustersOnTrack - 2))
        break;

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

          // Calculate the new and existing chi2 values
          double newchi2 =
              (bestTrack->chi2ndofZS() * bestTrack->chi2ndofZS() + bestTrack->chi2ndof() * bestTrack->chi2ndof());
          double oldchi2 = (conformalTrack->chi2ndofZS() * conformalTrack->chi2ndofZS() +
                            conformalTrack->chi2ndof() * conformalTrack->chi2ndof());

          //double deltachi2ZS = (bestTracks[itTrack]->chi2ndofZS() - conformalTracks[existingTrack]->chi2ndofZS());
          //double deltachi2   = (bestTracks[itTrack]->chi2ndof() - conformalTracks[existingTrack]->chi2ndof());

          // If the new track is an existing track + segment, make sure the increase in chi2 is not too much
          if (nOverlappingHits == conformalTrack->m_clusters.size()) {
            // Increase in chi2 is too much (double)
            if ((newchi2 - oldchi2) > oldchi2)
              break;

            // Otherwise replace the existing track
            conformalTrack = std::move(bestTrack);
            bestTrackUsed  = true;
          }

          // If the new track is a subtrack of an existing track, don't consider it further (already try removing bad hits from tracks
          else if (nOverlappingHits == bestTrack->m_clusters.size())
            break;

          // Otherwise take the longest if the delta chi2 is not too much
          else if (bestTrack->m_clusters.size() >= conformalTrack->m_clusters.size()) {  // New track longer/equal in length

            // Increase in chi2 is too much (double)
            if ((newchi2 - oldchi2) > oldchi2)
              break;

            // Otherwise take it
            conformalTrack = std::move(bestTrack);
            bestTrackUsed  = true;
          } else if (bestTrack->m_clusters.size() < conformalTrack->m_clusters.size()) {  // Old track longer

            // Must improve chi2 by factor two
            if ((newchi2 - 0.5 * oldchi2) > 0.)
              break;

            // Otherwise take it
            conformalTrack = std::move(bestTrack);
            bestTrackUsed  = true;
          }

          break;
        }
      }

      // If not a clone, save the new track
      if (!clone) {
        bestTrackUsed = true;
        if (debugSeed && kdhit == debugSeed) {
          streamlog_out(DEBUG7) << "== Pushing back best track with chi2/ndof " << bestTrack->chi2ndof() << std::endl;
        } else {
          streamlog_out(DEBUG7) << "Pushing back best track with chi2/ndof " << bestTrack->chi2ndof() << std::endl;
        }
        conformalTracks.push_back(std::move(bestTrack));
      }

      if (not bestTrackUsed) {
        bestTrack.reset();
      }

    }  // end for besttracks
  }
}

// Take a collection of tracks and try to extend them into the collection of clusters passed.
void ConformalTracking::extendTracks(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                     UKDTree& nearestNeighbours, Parameters const& parameters) {
  // Loop over all current tracks. At the moment this is a "stupid" algorithm: it will simply try to add every
  // hit in the collection to every track, and keep the ones thta have a good chi2. In fact, it will extrapolate
  // the track and do a nearest neighbours search, but this seemed to fail for some reason, TODO!

  streamlog_out(DEBUG7) << "EXTENDING " << conformalTracks.size() << " tracks, into " << collection.size() << " hits"
                        << std::endl;
  if (collection.size() == 0)
    return;
  int nTracks = conformalTracks.size();

  // Sort the hit collection by layer
  std::sort(collection.begin(), collection.end(), sort_by_layer);

  // First extend high pt tracks
  extendHighPT(conformalTracks, collection, nearestNeighbours, parameters);

  // Mark hits from "good" tracks as being used
  for (auto& conformalTrack : conformalTracks) {
    for (auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
  }

  for (int currentTrack = 0; currentTrack < nTracks; currentTrack++) {
    // Make sure that track hits are ordered from largest to smallest radius
    std::sort(conformalTracks[currentTrack]->m_clusters.begin(), conformalTracks[currentTrack]->m_clusters.end(),
              sort_by_radiusKD);

    // TODO: are kalman tracks necessary?
    // Make sure that all tracks have a kalman track attached to them
    //    if(conformalTracks[currentTrack]->kalmanTrack() == NULL){
    //      KalmanTrack* kalmanTrack = new KalmanTrack(conformalTracks[currentTrack]);
    //      conformalTracks[currentTrack]->setKalmanTrack(kalmanTrack);
    //    }

    // Get the associated MC particle
    MCParticle* associatedParticle = NULL;
    if (m_debugPlots) {
      associatedParticle = m_debugger.getAssociatedParticle(conformalTracks[currentTrack]);
      if (associatedParticle == NULL)
        streamlog_out(DEBUG7) << "- NULL particle!" << std::endl;

      // Check the track pt estimate
      TLorentzVector mc_helper;
      mc_helper.SetPxPyPzE(associatedParticle->getMomentum()[0], associatedParticle->getMomentum()[1],
                           associatedParticle->getMomentum()[2], associatedParticle->getEnergy());

      streamlog_out(DEBUG7) << "- extending track " << currentTrack << " with pt = " << mc_helper.Pt()
                            << ". pt estimate: " << conformalTracks[currentTrack]->pt() << " chi2/ndof "
                            << conformalTracks[currentTrack]->chi2ndof() << " and chi2/ndof ZS "
                            << conformalTracks[currentTrack]->chi2ndofZS() << std::endl;
    }

    // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
    int  nclusters = conformalTracks[currentTrack]->m_clusters.size();
    auto seedCell  = std::make_shared<Cell>(conformalTracks[currentTrack]->m_clusters[nclusters - 2],
                                           conformalTracks[currentTrack]->m_clusters[nclusters - 1]);

    // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
    //double searchDistance = parameters._maxDistance;
    //KDCluster* fakeHit = extrapolateCell(seedCell, searchDistance / 2.);  // TODO: make this search a function of radius
    //SharedKDClusters results2;
    //nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25 * searchDistance / 2., results2); //TEMP while not using
    //std::sort(results2.begin(), results2.end(), sort_by_radiusKD);
    //delete fakeHit;

    // Loop over all hits found and check if any have a sufficiently good delta Chi2
    SharedKDClusters goodHits;
    SKDCluster       bestCluster = NULL;
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

      // If this hit is on a new layer, then add the hit from the previous layer and start anew
      if (bestCluster != NULL && !(kdhit->sameLayer(bestCluster))) {
        bestCluster->used(true);
        conformalTracks[currentTrack]->add(bestCluster);
        conformalTracks[currentTrack]->linearRegression();
        conformalTracks[currentTrack]->linearRegressionConformal();
        bestCluster = NULL;
      }

      // Don't reuse hits
      if (kdhit->used()) {
        if (associated)
          streamlog_out(DEBUG7) << "used" << std::endl;
        continue;
      }

      // Don't pick up hits in the opposite side of the detector
      if ((conformalTracks[currentTrack]->m_clusters[nclusters - 1]->getZ() > 0. && kdhit->getZ() < 0.) ||
          (conformalTracks[currentTrack]->m_clusters[nclusters - 1]->getZ() < 0. && kdhit->getZ() > 0.)) {
        if (associated)
          streamlog_out(DEBUG7) << "opposite side of detector" << std::endl;
        continue;
      }

      // First check that the hit is not wildly away from the track (make cell and check angle)
      // SCell extensionCell = std::make_shared<Cell>(conformalTracks[currentTrack]->m_clusters[0],results2[newHit]);
      Cell   extensionCell(conformalTracks[currentTrack]->m_clusters[nclusters - 1], kdhit);
      double cellAngle   = seedCell->getAngle(extensionCell);
      double cellAngleRZ = seedCell->getAngleRZ(extensionCell);
      if (cellAngle > 3. * parameters._maxCellAngle || cellAngleRZ > 3. * parameters._maxCellAngleRZ) {
        if (associated)
          streamlog_out(DEBUG7) << "-- killed by cell angle cut" << std::endl;
        continue;
      }

      // Now fit the track with the new hit and check the increase in chi2
      double deltaChi2(0.), deltaChi2zs(0.);

      fitWithPoint(*(conformalTracks[currentTrack]), kdhit, deltaChi2,
                   deltaChi2zs);  //conformalTracks[currentTrack]->deltaChi2(results2[newHit]);

      if (associated) {
        streamlog_out(DEBUG7) << "-- hit was fitted and has a delta chi2 of " << deltaChi2 << " and delta chi2zs of "
                              << deltaChi2zs << std::endl;
      }

      // We have an estimate of the pT here, could use it in the chi2 criteria
      double chi2cut = parameters._chi2cut;

      //if (conformalTracks[currentTrack]->pt() < 5.)
      //  chi2cut = 1000.;

      if (deltaChi2 > chi2cut || deltaChi2zs > chi2cut)
        continue;

      if (associated)
        streamlog_out(DEBUG7) << "-- valid candidate!" << std::endl;

      if (bestCluster == NULL || deltaChi2 < bestChi2) {
        bestCluster = kdhit;
        bestChi2    = deltaChi2;
      }

      //conformalTracks[currentTrack]->add(kdhit);
      //kdhit->used(true);

      // New hit has been added, go back to the beginning and start again...
      // nKDHit = 0;
    }

    if (bestCluster != NULL) {
      bestCluster->used(true);
      conformalTracks[currentTrack]->add(bestCluster);
      conformalTracks[currentTrack]->linearRegression();
      conformalTracks[currentTrack]->linearRegressionConformal();
      bestCluster = NULL;
    }

  }  // End of loop over tracks
}

// Extend seed cells
void ConformalTracking::extendSeedCells(SharedCells& cells, UKDTree& nearestNeighbours, bool extendingTrack,
                                        const SharedKDClusters& /*debugHits*/, Parameters const& parameters,
                                        bool vertexToTracker) {
  unsigned int nCells   = 0;
  int          depth    = 0;
  int          startPos = 0;

  // Keep track of existing cells in case there are branches in the track
  std::map<SKDCluster, SharedCells> existingCells;

  // Try to create all "downstream" cells until no more can be added
  while (cells.size() != nCells) {
    // Extend all cells with depth N. In the next iteration, look at cells with depth N+1
    nCells = cells.size();
    for (unsigned int itCell = startPos; itCell < nCells; itCell++) {
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
            if (hit->sameLayer(nhit))
              return true;
            if (nhit->endcap() && hit->endcap() && (nhit->forward() != hit->forward()))
              return true;
            if ((vertexToTracker && nhit->getR() >= hit->getR()) || (!vertexToTracker && nhit->getR() <= hit->getR()))
              return true;
            return false;
          });

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

        if (extendingTrack)
          streamlog_out(DEBUG7) << "looking at neighbour " << neighbour << " at u,v: " << nhit->getU() << "," << nhit->getV()
                                << std::endl;

        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
        //        if(used.count(nhit)){if(extendingTrack)streamlog_out(DEBUG7)<<"- used"<<std::endl; continue;}
        if (nhit->used())
          continue;
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
              alreadyExists = true;

              // Check if cell angle is too large to rejoin
              if (cells[itCell]->getAngle(existingCell) > parameters._maxCellAngle ||
                  cells[itCell]->getAngleRZ(existingCell) > parameters._maxCellAngleRZ)
                continue;

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

  // No more downstream cells can be added
}

void ConformalTracking::extendHighPT(UniqueKDTracks& conformalTracks, SharedKDClusters& /*collection*/,
                                     UKDTree& nearestNeighbours, Parameters const& parameters, bool /*radialSearch*/) {
  // Set the cell angle criteria for high pt tracks
  Parameters highPTParameters(parameters);
  highPTParameters._maxCellAngle   = 0.005;
  highPTParameters._maxCellAngleRZ = 0.005;

  // Loop over all tracks
  int nTracks = conformalTracks.size();
  for (int currentTrack = 0; currentTrack < nTracks; currentTrack++) {
    // Get the track
    UKDTrack& track = conformalTracks[currentTrack];

    // Only look at high pt tracks
    if (track->pt() < 10.)
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
    UniqueCellularTracks extensions;

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

      UKDTrack extendedTrack = std::unique_ptr<KDTrack>(new KDTrack());
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
      conformalTracks[currentTrack] = std::move(lowestChi2Track);
    }
    streamlog_out(DEBUG7) << "- track finishes with size " << conformalTracks[currentTrack]->m_clusters.size() << std::endl;
  }

  return;
}

//===================================
// Cellular Track Functions
//===================================

// New test at creating cellular tracks. In this variant, don't worry about clones etc, give all possible routes back to the seed cell. Then cut
// on number of clusters on each track, and pass back (good tracks to then be decided based on best chi2
void ConformalTracking::createTracksNew(UniqueCellularTracks& finalcellularTracks, SCell& seedCell,
                                        std::map<SCell, bool>& /*usedCells*/) {
  // Final container to be returned
  UniqueCellularTracks cellularTracks;

  // Make the first cellular track using the seed cell
  auto seedTrack = std::unique_ptr<cellularTrack>(new cellularTrack());
  seedTrack->push_back(seedCell);
  cellularTracks.push_back(std::move(seedTrack));

  // Now start to follow all paths back from this seed cell
  // While there are still tracks that are not finished (last cell weight 0), keep following their path
  while (toBeUpdated(cellularTracks)) {
    //   streamlog_out(DEBUG7)<<"== Updating "<<cellularTracks.size()<<" tracks"<<std::endl;
    // Loop over all (currently existing) tracks
    int nTracks = cellularTracks.size();
    if (nTracks > 10000) {
      streamlog_out(WARNING) << "Going to create " << nTracks << std::endl;
    }
    if (nTracks > 5e5) {
      streamlog_out(ERROR) << "Too many tracks (" << nTracks << " > 1e6) are going to be created, tightening parameters"
                           << std::endl;
      throw TooManyTracksException(this);
    }
    for (int itTrack = 0; itTrack < nTracks; itTrack++) {
      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom().size() == 0) {
        //       streamlog_out(DEBUG7)<<"-- Track "<<itTrack<<" is finished"<<std::endl;
        continue;
      }

      // While there is only one path leading from this cell, follow that path
      SCell cell = cellularTracks[itTrack]->back();
      //     streamlog_out(DEBUG7)<<"-- Track "<<itTrack<<" has "<<(*(cell->getFrom())).size()<<" cells attached to the end of it"<<std::endl;
      //      while(cell->getWeight() > 0 && (*(cell->getFrom())).size() == 1){
      while (cell->getFrom().size() == 1) {
        //       streamlog_out(DEBUG7)<<"- simple extension"<<std::endl;
        // Get the cell that it attaches to
        auto parentCell = SCell(cell->getFrom()[0]);
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

      //     streamlog_out(DEBUG7)<<"- making "<<nBranches<<" branches"<<std::endl;

      // For each additional branch make a new track
      for (int itBranch = 1; itBranch < nBranches; itBranch++) {
        auto branchedTrack = std::unique_ptr<cellularTrack>(new cellularTrack(*(cellularTracks[itTrack].get())));
        branchedTrack->push_back(SCell(cell->getFrom()[itBranch]));
        cellularTracks.push_back(std::move(branchedTrack));
      }

      // Keep the existing track for the first branch
      cellularTracks[itTrack]->push_back(SCell(cell->getFrom()[0]));
    }
  }

  finalcellularTracks.insert(finalcellularTracks.end(), std::make_move_iterator(cellularTracks.begin()),
                             std::make_move_iterator(cellularTracks.end()));

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
void ConformalTracking::getFittedTracks(UniqueKDTracks& finalTracks, UniqueCellularTracks& candidateTracks,
                                        std::map<SCell, bool>& /*usedCells*/, Parameters const& parameters) {
  // Make a container for all tracks being considered, initialise variables
  UniqueKDTracks trackContainer;
  //  std::vector<double> trackChi2ndofs;

  // Loop over all candidate tracks and do an inital fit to get the track angle (needed to calculate the
  // hit errors for the error-weighted fit)
  for (auto& candidateTrack : candidateTracks) {
    // If there are not enough hits on the track, ignore it
    if (int(candidateTrack->size()) < (parameters._minClustersOnTrack - 2)) {
      candidateTrack.reset();
      continue;
    }

    // Make the fitting object. TGraphErrors used for 2D error-weighted fitting
    UKDTrack track = std::unique_ptr<KDTrack>(new KDTrack());

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
    track->linearRegressionConformal();               //FCC study
    double chi2ndof = track->chi2() / (npoints - 2);  //FCC study

    // We try to see if there are spurious hits causing the chi2 to be very large. This would cause us to throw away
    // good tracks with perhaps just a single bad hit. Try to remove each hit and see if fitting without it causes a
    // significant improvement in chi2/ndof

    // Loop over each hit (starting at the back, since we will use the 'erase' function to get rid of them)
    // and see if removing it improves the chi2/ndof
    double removed = 0;
    if (chi2ndof > parameters._chi2cut &&
        chi2ndof <
            parameters
                ._chi2cut) {  //CHANGE ME?? Upper limit to get rid of really terrible tracks (temp lower changed from 0 to parameters._chi2cut)
      for (int point = npoints - 1; point >= 0; point--) {
        // Stop if we would remove too many points on the track to meet the minimum hit requirement (or if the track has more
        // than 2 hits removed)
        if ((npoints - removed - 1) < parameters._minClustersOnTrack || removed == 2)
          break;

        // Refit the track without this point
        double newChi2ndof = fitWithoutPoint(*track, point);

        // If the chi2/ndof is significantly better, remove the point permanently CHANGE ME??
        //        if( (chi2ndof - newChi2ndof) > 0 && (chi2ndof - newChi2ndof) > 1. ){
        if ((newChi2ndof - chi2ndof) < chi2ndof) {
          track->remove(point);
          removed++;
          chi2ndof = newChi2ndof;
        }
      }
    }

    // Store the track information
    trackContainer.push_back(std::move(track));
    //    trackChi2ndofs.push_back(chi2ndof);

    candidateTrack.reset();

  }  // end for candidateTracks

  // Now have all sets of conformal tracks and their chi2/ndof. Decide which tracks to send back, ie. the one with
  // lowest chi2/ndof, and possibly others if they are not clones and have similar chi2 value
  getLowestChi2(finalTracks, trackContainer);

  // Send back the final set of tracks
  return;
}

// Pick the lowest chi2/ndof KDTrack from a list of possible tracks, and additionally return other tracks in the collection with similar chi2/ndof values that don't share many hits
void ConformalTracking::getLowestChi2(UniqueKDTracks& finalTracks, UniqueKDTracks& trackContainer) {
  // Get the lowest chi2/ndof value from the given tracks
  //  double lowestChi2ndof = *std::min_element(trackChi2ndofs.begin(),trackChi2ndofs.end());
  UKDTrack& lowestChi2ndofTrack = trackContainer[0];
  double    lowestChi2ndof      = sqrt(lowestChi2ndofTrack->chi2ndof() * lowestChi2ndofTrack->chi2ndof() +
                               lowestChi2ndofTrack->chi2ndofZS() * lowestChi2ndofTrack->chi2ndofZS());

  for (auto& itTrack : trackContainer) {
    double chi2ndof = sqrt(itTrack->chi2ndof() * itTrack->chi2ndof() + itTrack->chi2ndofZS() * itTrack->chi2ndofZS());
    if (chi2ndof < lowestChi2ndof) {
      lowestChi2ndof = itTrack->chi2ndof();
      // lowestChi2ndofTrack = trackContainer[itTrack];
    }
  }

  //finalTracks.push_back(lowestChi2ndofTrack);
  //return;

  // Loop over all other tracks and decide whether or not to save them
  // Look at the difference in chi2/ndof - we want to keep tracks with similar chi2/ndof. If they
  // are clones then take the longest
  copy_if(std::make_move_iterator(trackContainer.begin()), std::make_move_iterator(trackContainer.end()),
          std::back_inserter(finalTracks), [lowestChi2ndof](UKDTrack const& track) {
            return ((sqrt((track->chi2ndof() * track->chi2ndof() + track->chi2ndofZS() * track->chi2ndofZS())) -
                     lowestChi2ndof) < 10.);
          });
  trackContainer.clear();

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
  //int npoints = track.nPoints();

  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal();

  double chi2ndof = sqrt(track.chi2ndofZS() * track.chi2ndofZS() + track.chi2ndof() * track.chi2ndof());
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
  track.linearRegressionConformal();  // FCC study

  double chi2ndof   = track.chi2ndof();
  double chi2ndofZS = track.chi2ndofZS();  // FCC study

  return sqrt(chi2ndof * chi2ndof + chi2ndofZS * chi2ndofZS);
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
double ConformalTracking::checkReal(UKDTrack& track, std::map<SKDCluster, MCParticle*> kdParticles,
                                    std::map<MCParticle*, bool>&            reconstructed,
                                    std::map<MCParticle*, SharedKDClusters> MCparticleHits) {
  // Store all mcparticles associated to this track
  std::vector<MCParticle*> particles;
  std::map<MCParticle*, double> particleHits;
  double nHits = 0.;

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
  MCParticle* bestParticle = NULL;
  for (auto& particle : particles) {
    if (particleHits[particle] > bestHits) {
      bestHits     = particleHits[particle];
      bestParticle = particle;
    }
  }

  // Calculate the purity
  double purity = bestHits / nHits;
  streamlog_out(DEBUG7) << "Number of hits on track: " << nHits << ". Good hits: " << bestHits << ". Purity: " << purity
                        << ". Pt: " << sqrt(bestParticle->getMomentum()[0] * bestParticle->getMomentum()[0] +
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
    streamlog_out(DEBUG7) << "Hit " << th << " u = " << trackHits[th]->getU() << " v = " << trackHits[th]->getV()
                          << " x = " << trackHits[th]->getX() << " y = " << trackHits[th]->getY()
                          << " z = " << trackHits[th]->getZ() << std::endl;

    // Check which particle the hit belongs to
    MCParticle* particle  = kdParticles[trackHits[th]];
    double      mcVertexX = particle->getVertex()[0];
    double      mcVertexY = particle->getVertex()[1];
    double      mcVertexR = sqrt(pow(mcVertexX, 2) + pow(mcVertexY, 2));

    streamlog_out(DEBUG7) << "Come from particle with id " << particle->getPDG() << " produced at vertexR = " << mcVertexR;

    if (particle->getParents().size() != 0) {
      streamlog_out(DEBUG7) << ". Comes from particle with id " << particle->getParents()[0]->getPDG();
    }
    streamlog_out(DEBUG7) << std::endl;
  }

  //  for(int itCluster=0;itCluster<clusters.size();itCluster++)streamlog_out(DEBUG7)<<"Hit "<<itCluster<<" has position: "<<clusters[itCluster]->getU()<<","<<clusters[itCluster]->getV()<<std::endl;
  streamlog_out(DEBUG7) << "== Terms in conformal fit: " << track->m_intercept << ", " << track->m_gradient << ", "
                        << track->m_quadratic << std::endl;
  streamlog_out(DEBUG7) << "== Terms in zs fit: " << track->m_interceptZS << ", " << track->m_gradientZS << std::endl;

  track->linearRegressionConformal(true);
  track->calculateChi2SZ(NULL, true);
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
void ConformalTracking::checkReconstructionFailure(MCParticle* particle,
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
  UniqueCellularTracks cellularTracks;
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
  mcTrack->calculateChi2SZ(NULL, true);

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
