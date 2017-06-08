#include "ConformalTracking.h"
#include "DDRec/API/IDDecoder.h"

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

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/LCDD.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

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

#include "Cell.h"
#include "KDTree.h"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace DD4hep;
using namespace AIDA;

// Static instantiation of the processor
ConformalTracking aConformalTracking;

/*
 
 Pattern recognition code for the CLIC detector, using conformal mapping and cellular automaton
 
 */

ConformalTracking::ConformalTracking() : Processor("ConformalTracking") {
  // Processor description
  _description = "ConformalTracking constructs tracks using a combined conformal mapping and cellular automaton approach.";

  // Input collections - tracker hits
  std::vector<std::string> inputTrackerHitCollections;
  inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("VXDEndcapTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerEndcapHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerEndcapHits"));
  registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames", "Name of the TrackerHit input collections",
                           m_inputTrackerHitCollections, inputTrackerHitCollections);

  // Debugging collections - MC particles and relation collections
  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollectionName", "Name of the MCParticle input collection",
                          m_inputParticleCollection, std::string("MCParticle"));
  std::vector<std::string> inputRelationCollections;
  inputRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
  inputRelationCollections.push_back(std::string("VXDEndcapTrackerHitRelations"));
  inputRelationCollections.push_back(std::string("InnerTrackerBarrelHitsRelations"));
  inputRelationCollections.push_back(std::string("OuterTrackerBarrelHitsRelations"));
  inputRelationCollections.push_back(std::string("InnerTrackerEndcapHitsRelations"));
  inputRelationCollections.push_back(std::string("OuterTrackerEndcapHitsRelations"));
  registerInputCollections(LCIO::LCRELATION, "RelationsNames", "Name of the TrackerHit relation collections",
                           m_inputRelationCollections, inputRelationCollections);

  // Output collections - tracks
  registerOutputCollection(LCIO::TRACK, "SiTrackCollectionName", "Silicon track Collection Name", m_outputTrackCollection,
                           std::string("CATracks"));
  registerOutputCollection(LCIO::TRACKERHITPLANE, "DebugHits", "DebugHits", m_outputDebugHits, std::string("DebugHits"));

  // Parameters for tracking
  registerProcessorParameter("DebugPlots", "Plots for debugging the tracking", m_debugPlots, bool(false));
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
  DD4hep::Geometry::LCDD& lcdd        = DD4hep::Geometry::LCDD::getInstance();
  const double            position[3] = {0, 0, 0};  // position to calculate magnetic field at (the origin in this case)
  double                  magneticFieldVector[3] = {0, 0, 0};  // initialise object to hold magnetic field
  lcdd.field().magneticField(position, magneticFieldVector);   // get the magnetic field vector from DD4hep
  m_magneticField = magneticFieldVector[2] / dd4hep::tesla;    // z component at (0,0,0)

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

  streamlog_out(DEBUG4) << "Event number: " << m_eventNumber << std::endl;

  // Set up ID decoder
  UTIL::BitField64 m_encoder(lcio::LCTrackerCellID::encoding_string());

  // Object to store all of the track hit collections passed to the pattern recognition
  std::vector<LCCollection*>        trackerHitCollections;
  std::vector<LCRelationNavigator*> relations;
  LCCollection*                     particleCollection;

  // Loop over each input collection and get the hits
  for (unsigned int collection = 0; collection < m_inputTrackerHitCollections.size(); collection++) {
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt);
    if (trackerHitCollection == 0)
      continue;
    streamlog_out(DEBUG4) << "Collection " << m_inputTrackerHitCollections[collection] << " contains "
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
      LCRelationNavigator* relation = new LCRelationNavigator(trackerHitRelationCollection);
      relations.push_back(relation);
    }
  }

  // Get the MC particle collection
  getCollection(particleCollection, m_inputParticleCollection, evt);

  // Make the output track collection
  LCCollectionVec* trackCollection    = new LCCollectionVec(LCIO::TRACK);
  LCCollectionVec* debugHitCollection = new LCCollectionVec(LCIO::TRACKERHITPLANE);
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
  std::map<int, std::vector<KDCluster*>> collectionClusters;  // Conformal hits
  std::map<KDCluster*, TrackerHitPlane*> kdClusterMap;        // Their link to "real" hits
  std::map<TrackerHitPlane*, KDCluster*> conformalHits;       // The reverse link

  // Debug collections (not filled if debug off)
  std::map<KDCluster*, MCParticle*>              kdParticles;    // Link from conformal hit to MC particle
  std::map<MCParticle*, std::vector<KDCluster*>> particleHits;   // List of conformal hits on each MC particle
  std::map<MCParticle*, bool>                    reconstructed;  // Check for MC particles
  std::vector<KDCluster*> debugHits;                             // Debug hits for plotting

  // Create the conformal hit collections for each tracker hit collection (and save the link)
  for (unsigned int collection = 0; collection < trackerHitCollections.size(); collection++) {
    // Loop over tracker hits and make conformal hit collection
    std::vector<KDCluster*> tempClusters;
    int                     nHits = trackerHitCollections[collection]->getNumberOfElements();
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
      if (side != ILDDetID::barrel)
        isEndcap = true;

      // Make a new kd cluster
      KDCluster* kdhit = new KDCluster(hit, isEndcap);

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
        // Store the information
        particleHits[particle].push_back(kdhit);
        kdParticles[kdhit] = particle;
        // Draw plots for event 0
        if (m_eventNumber == 0) {
          m_conformalEvents->Fill(kdhit->getU(), kdhit->getV());
          m_nonconformalEvents->Fill(hit->getPosition()[0], hit->getPosition()[1]);
          m_conformalEventsRTheta->Fill(kdhit->getR(), kdhit->getTheta());
        }
      }
    }
    collectionClusters[collection] = tempClusters;
  }

  // WHAT TO DO ABOUT THIS?? POSSIBLY MOVE DEPENDING ON MC RECONSTRUCTION (and in fact, would fit better into the check reconstruction code at present)

  // Now loop over all MC particles and make the cells connecting hits
  int nParticles = particleCollection->getNumberOfElements();
  for (int itP = 0; itP < nParticles; itP++) {
    // Get the particle
    MCParticle* mcParticle = dynamic_cast<MCParticle*>(particleCollection->getElementAt(itP));
    // Get the vector of hits from the container
    if (particleHits.count(mcParticle) == 0)
      continue;
    std::vector<KDCluster*> trackHits = particleHits[mcParticle];
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
    KDTrack* mcTrack = new KDTrack();
    // Loop over all hits for debugging
    for (int itHit = 0; itHit < trackHits.size(); itHit++) {
      // Get the conformal clusters
      KDCluster* cluster = trackHits[itHit];
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
    delete mcTrack;

    // Now loop over the hits and make cells - filling histograms along the way
    int nHits = trackHits.size();
    for (int itHit = 0; itHit < (nHits - 2); itHit++) {
      // Get the conformal clusters
      KDCluster* cluster0 = trackHits[itHit];
      KDCluster* cluster1 = trackHits[itHit + 1];
      KDCluster* cluster2 = trackHits[itHit + 2];

      // Make the two cells connecting these three hits
      Cell* cell = new Cell(cluster0, cluster1);
      cell->setWeight(itHit);
      Cell* cell1 = new Cell(cluster1, cluster2);
      cell1->setWeight(itHit + 1);

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
      delete cell;
      delete cell1;
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
  }

  // END OF "WHAT TO DO ABOUT THIS??"

  //------------------------------------------------------------------------------
  // Now the track reconstruction strategy. Perform a sequential search, with hits
  // removed from the seeding collections once tracks have been built
  //------------------------------------------------------------------------------

  // The final vector of conformal tracks
  std::vector<KDTrack*>   conformalTracks;
  std::vector<KDCluster*> kdClusters;
  KDTree*                 nearestNeighbours = NULL;

  // Build tracks in the vertex barrel
  std::vector<int> vertexHits = {0};
  combineCollections(kdClusters, nearestNeighbours, vertexHits, collectionClusters);
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours);

  // Mark hits from "good" tracks as being used
  for (unsigned int itTrack = 0; itTrack < conformalTracks.size(); itTrack++) {
    for (unsigned int itHit = 0; itHit < conformalTracks[itTrack]->m_clusters.size(); itHit++)
      conformalTracks[itTrack]->m_clusters[itHit]->used(true);
  }

  // Extend through the endcap
  std::vector<int> vertexEndcapHits = {1};
  combineCollections(kdClusters, nearestNeighbours, vertexEndcapHits, collectionClusters);
  extendTracks(conformalTracks, kdClusters, nearestNeighbours);

  // Make combined vertex tracks
  std::vector<int> vertexCombinedHits = {0, 1};
  combineCollections(kdClusters, nearestNeighbours, vertexCombinedHits, collectionClusters);
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours);

  // Mark hits from "good" tracks as being used
  for (unsigned int itTrack = 0; itTrack < conformalTracks.size(); itTrack++) {
    for (unsigned int itHit = 0; itHit < conformalTracks[itTrack]->m_clusters.size(); itHit++)
      conformalTracks[itTrack]->m_clusters[itHit]->used(true);
  }

  // Make leftover tracks in the vertex with lower requirements

  // Store global variables
  double maxCellAngle       = m_maxCellAngle;
  double maxCellAngleRZ     = m_maxCellAngleRZ;
  double chi2cut            = m_chi2cut;
  double minClustersOnTrack = m_minClustersOnTrack;

  // Lower cell angle criteria
  m_maxCellAngle *= 20.;
  m_maxCellAngleRZ *= 20.;
  //  m_chi2cut*=20.;
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours);

  // Mark hits from "good" tracks as being used
  for (unsigned int itTrack = 0; itTrack < conformalTracks.size(); itTrack++) {
    for (unsigned int itHit = 0; itHit < conformalTracks[itTrack]->m_clusters.size(); itHit++)
      conformalTracks[itTrack]->m_clusters[itHit]->used(true);
  }

  // Lower number of hits on track
  m_minClustersOnTrack = 4;
  buildNewTracks(conformalTracks, kdClusters, nearestNeighbours);

  // Mark hits from "good" tracks as being used
  for (unsigned int itTrack = 0; itTrack < conformalTracks.size(); itTrack++) {
    for (unsigned int itHit = 0; itHit < conformalTracks[itTrack]->m_clusters.size(); itHit++)
      conformalTracks[itTrack]->m_clusters[itHit]->used(true);
  }

  // Put back original parameters
  m_maxCellAngle       = maxCellAngle;
  m_maxCellAngleRZ     = maxCellAngleRZ;
  m_chi2cut            = chi2cut;
  m_minClustersOnTrack = minClustersOnTrack;

  // Sort by pt (low to hight)
  std::sort(conformalTracks.begin(), conformalTracks.end(), sort_by_pt);

  /*
  // Extend them through the inner and outer trackers
  std::vector<int> trackerHits = {2, 3, 4, 5};
  combineCollections(kdClusters, nearestNeighbours, trackerHits, collectionClusters);
  extendTracks(conformalTracks, kdClusters, nearestNeighbours);

  // Increase chi2 to add hits
  double chi2increase = m_chi2increase;
  m_chi2increase      = 1000.;
  extendTracks(conformalTracks, kdClusters, nearestNeighbours);
  m_chi2increase = chi2increase;
  */

  // Clean up
  delete nearestNeighbours;

  // Now in principle have all conformal tracks, but due to how the check for clones is performed (ish) there is a possibility
  // that clones/fakes are still present. Try to remove them by looking at overlapping hits. Turned off at the moment

  std::vector<KDTrack*> conformalTracksFinal = conformalTracks;

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
  streamlog_out(DEBUG) << "*** CA has made " << conformalTracksFinal.size()
                       << (conformalTracksFinal.size() == 1 ? " track ***" : " tracks ***") << std::endl;

  // Loop over all track candidates
  for (unsigned int caTrack = 0; caTrack < conformalTracksFinal.size(); caTrack++) {
    // Vector of all the hits on the track
    KDTrack* conformalTrack = conformalTracksFinal[caTrack];
    streamlog_out(DEBUG5) << "Made a track with " << conformalTrack->m_clusters.size() << " hits" << std::endl;

    // Make the LCIO track hit vector
    EVENT::TrackerHitVec trackHits;
    for (unsigned int itHit = 0; itHit < conformalTrack->m_clusters.size(); itHit++) {
      KDCluster* cluster = conformalTrack->m_clusters[itHit];
      trackHits.push_back(kdClusterMap[cluster]);
    }
    // Add kalman filtered hits
    if (conformalTrack->kalmanTrack() != NULL) {
      KalmanTrack* kalmanTrack = conformalTrack->kalmanTrack();
      for (int i = 1; i < kalmanTrack->m_kalmanClusters.size(); i++) {
        KDCluster* cluster = kalmanTrack->m_kalmanClusters[i];
        trackHits.push_back(kdClusterMap[cluster]);
      }
    }

    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(), trackHits.end(), sort_by_radius);

    // Now we can make the track object and relations object, and fit the track
    TrackImpl* track = new TrackImpl;

    // First, for some reason there are 2 track objects, one which gets saved and one which is used for fitting. Don't ask...
    MarlinTrk::IMarlinTrack* marlinTrack = trackFactory->createTrack();

    // Make an initial covariance matrix with very broad default values
    EVENT::FloatVec covMatrix(15, 0);             // Size 15, filled with 0s
    covMatrix[0]  = (m_initialTrackError_d0);     //sigma_d0^2
    covMatrix[2]  = (m_initialTrackError_phi0);   //sigma_phi0^2
    covMatrix[5]  = (m_initialTrackError_omega);  //sigma_omega^2
    covMatrix[9]  = (m_initialTrackError_z0);     //sigma_z0^2
    covMatrix[14] = (m_initialTrackError_tanL);   //sigma_tanl^2

    // Try to fit
    int fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackHits, track, MarlinTrk::IMarlinTrack::forward,
                                                       covMatrix, m_magneticField, m_maxChi2perHit);

    // Check track quality - if fit fails chi2 will be 0. For the moment add hits by hand to any track that fails the track fit, and store it as if it were ok...
    if (track->getChi2() == 0.) {
      streamlog_out(DEBUG) << "Fit failed. Track has " << track->getTrackerHits().size() << " hits" << std::endl;
      streamlog_out(DEBUG) << "Fit fail error " << fitError << std::endl;
      // for(unsigned int p=0;p<trackHits.size();p++){
      //   track->addHit(trackHits[p]);
      // }
      //    }//
      //    delete marlinTrack;
      delete track;
      delete marlinTrack;
      continue;
    }

    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2] = trackHits.size();

    // Push back to the output container
    trackCollection->addElement(track);
  }

  // calculate purities and check if tracks have been reconstructed
  if (m_debugPlots) {
    for (int itrack = 0; itrack < conformalTracksFinal.size(); itrack++) {
      KDTrack* debugTrack = conformalTracksFinal[itrack];

      m_conformalChi2->Fill(debugTrack->chi2ndof());
      streamlog_out(DEBUG) << "-------------------- New TRACK --------------------" << std::endl;
      double purity = checkReal(debugTrack, kdParticles, reconstructed, particleHits);
      if (purity >= m_purity) {
        m_conformalChi2real->Fill(debugTrack->chi2ndof());
      }
      if (purity < m_purity) {
        m_conformalChi2fake->Fill(debugTrack->chi2ndof());
      }
      m_conformalChi2Purity->Fill(purity, debugTrack->chi2ndof());
    }
  }  //*/

  // Draw the cells for all produced tracks
  if (m_debugPlots && m_eventNumber == 0) {
    m_canvConformalEventDisplay->cd();
    for (int itrack = 0; itrack < conformalTracksFinal.size(); itrack++) {
      KDTrack*                debugTrack = conformalTracksFinal[itrack];
      std::vector<KDCluster*> clusters   = debugTrack->m_clusters;
      for (int itCluster = 1; itCluster < clusters.size(); itCluster++)
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
    vector<KDCluster*> kdClusters;
    for (int i = 0; i < trackerHitCollections.size(); i++)
      kdClusters.insert(kdClusters.begin(), collectionClusters[i].begin(), collectionClusters[i].end());
    KDTree* nearestNeighbours = new KDTree(kdClusters, m_thetaRange);

    for (int itP = 0; itP < nParticles; itP++) {
      // Get the particle
      MCParticle* mcParticle = dynamic_cast<MCParticle*>(particleCollection->getElementAt(itP));
      // Get the conformal hits
      if (particleHits.count(mcParticle) == 0)
        continue;
      std::vector<KDCluster*> mcHits = particleHits[mcParticle];
      // Cut on the number of hits
      int uniqueHits = getUniqueHits(mcHits);
      if (uniqueHits < m_minClustersOnTrack)
        continue;
      // Check if it was stable
      if (mcParticle->getGeneratorStatus() != 1)
        continue;
      // Check if it was reconstructed
      streamlog_out(DEBUG) << "-------------------- New PARTICLE --------------------" << std::endl;
      // List the pt
      double particlePt = sqrt(mcParticle->getMomentum()[0] * mcParticle->getMomentum()[0] +
                               mcParticle->getMomentum()[1] * mcParticle->getMomentum()[1]);
      streamlog_out(DEBUG) << "Particle pt: " << particlePt << std::endl;
      checkReconstructionFailure(mcParticle, particleHits, nearestNeighbours);
      if (reconstructed.count(mcParticle)) {
        nReconstructed++;
        continue;
      }
      // Draw the cells connecting the hits
      std::sort(mcHits.begin(), mcHits.end(), sort_by_radiusKD);
      for (int iHit = 0; iHit < (mcHits.size() - 1); iHit++) {
        drawline(mcHits[iHit], mcHits[iHit + 1], iHit + 1);
      }
      streamlog_out(DEBUG) << "Unreconstructed particle pt: " << particlePt << std::endl;
      nUnreconstructed++;

      // Check why particles were not reconstructed
      //      checkReconstructionFailure(mcParticle, particleHits, used, nearestNeighbours);
    }
    streamlog_out(DEBUG) << "Reconstructed " << nReconstructed << " particles out of " << nReconstructed + nUnreconstructed
                         << ". Gives efficiency "
                         << 100. * (double)nReconstructed / (double)(nReconstructed + nUnreconstructed) << "%" << std::endl;
    delete nearestNeighbours;
  }

  // Save the output track collection
  evt->addCollection(trackCollection, m_outputTrackCollection);
  evt->addCollection(debugHitCollection, m_outputDebugHits);

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
void ConformalTracking::combineCollections(std::vector<KDCluster*>& kdClusters, KDTree*& nearestNeighbours,
                                           std::vector<int> combination,
                                           std::map<int, std::vector<KDCluster*>> collectionClusters) {
  // Clear the input objects
  kdClusters.clear();
  if (nearestNeighbours != NULL)
    delete nearestNeighbours;

  // Loop over all given collections
  for (unsigned int i = 0; i < combination.size(); i++) {
    // Copy the clusters to the output vector
    std::vector<KDCluster*> clusters = collectionClusters[combination[i]];  //this makes a copy FIX ME
    int                     nhits    = clusters.size();
    for (int hit = 0; hit < nhits; hit++) {
      kdClusters.push_back(clusters[hit]);
    }
  }

  // Sort the KDClusters from larger to smaller radius
  std::sort(kdClusters.begin(), kdClusters.end(), sort_by_radiusKD);

  // Make the binary search tree. This tree class contains two binary trees - one sorted by u-v and the other by theta
  nearestNeighbours = new KDTree(kdClusters, m_thetaRange);
}

// Take a collection of hits and try to produce tracks out of them
void ConformalTracking::buildNewTracks(std::vector<KDTrack*>& conformalTracks, std::vector<KDCluster*>& collection,
                                       KDTree* nearestNeighbours, bool radialSearch) {
  streamlog_out(DEBUG) << "BUILDING new tracks" << std::endl;

  // Sort the input collection by radius
  std::sort(collection.begin(), collection.end(), sort_by_radiusKD);

  // Loop over all hits, using each as a seed to produce a new track
  unsigned int nKDHits = collection.size();
  for (unsigned int nKDHit = 0; nKDHit < nKDHits; nKDHit++) {
    // Get the kdHit and check if it has already been used (assigned to a track)
    KDCluster* kdhit = collection[nKDHit];
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG) << "Starting to seed with debug cluster" << std::endl;
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
    VecCluster results;
    double     theta = kdhit->getTheta();
    if (radialSearch)
      nearestNeighbours->allNeighboursInRadius(kdhit, m_maxDistance, results);
    else
      nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);

    // Sort the neighbours from outer to inner radius
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG) << "- picked up " << results.size() << " neighbours from theta search" << std::endl;
    if (results.size() == 0)
      continue;
    std::sort(results.begin(), results.end(), sort_by_radiusKD);

    // Objects to hold cells
    std::vector<Cell*> cells;

    // Make seed cells pointing inwards (conformal space)
    for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
      // Get the neighbouring hit
      KDCluster* nhit = results[neighbour];

      // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
      if (nhit->used())
        continue;
      if (kdhit->sameLayer(nhit))
        continue;
      if (nhit->getR() >= kdhit->getR())
        continue;

      // Check if the cell would be too long (hit very far away)
      double length = sqrt((kdhit->getU() - nhit->getU()) * (kdhit->getU() - nhit->getU()) +
                           (kdhit->getV() - nhit->getV()) * (kdhit->getV() - nhit->getV()));
      if (length > m_maxDistance)
        continue;

      // Create the new seed cell
      Cell* cell = new Cell(kdhit, nhit);
      cells.push_back(cell);
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG) << "- made cell with neighbour " << neighbour << " at " << nhit->getU() << "," << nhit->getV()
                             << std::endl;

      // Debug plotting
      if (m_debugPlots && m_eventNumber == 0) {
        m_canvConformalEventDisplayAllCells->cd();
        drawline(kdhit, nhit, 1);
      }
    }

    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG) << "- produced " << cells.size() << " seed cells" << std::endl;

    // No seed cells produced
    if (cells.size() == 0)
      continue;

    // All seed cells have been created, now try create all "downstream" cells until no more can be added
    std::vector<KDCluster*> debugHits;
    extendSeedCells(cells, nearestNeighbours, false, debugHits);

    // Now have all cells stemming from this seed hit. If it is possible to produce a track (ie. cells with depth X) then we will now...
    //      if(depth < (m_minClustersOnTrack-1)) continue; // TODO: check if this is correct

    // We create all acceptable tracks by looping over all cells with high enough weight to create
    // a track and trace their route back to the seed hit. We then have to choose the best candidate
    // at the end (by minimum chi2 of a linear fit)
    std::map<Cell*, bool> usedCells;
    std::map<Cell*, bool> usedCells2;
    std::vector<KDTrack*> cellTracks;

    // Sort Cells from highest to lowest weight
    std::sort(cells.begin(), cells.end(), sort_by_cellWeight);

    // Create tracks by following a path along cells
    int nCells = cells.size();
    for (int itCell = 0; itCell < nCells; itCell++) {
      // Check if this cell has already been used
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG) << "-- looking at cell " << itCell << std::endl;
      if (usedCells.count(cells[itCell]))
        continue;

      // Check if this cell could produce a track (is on a long enough chain)
      if (cells[itCell]->getWeight() < (m_minClustersOnTrack - 2))
        break;

      // Produce all tracks leading back to the seed hit from this cell
      std::vector<cellularTrack*> candidateTracks;
      std::vector<cellularTrack*> candidateTracksTemp;
      createTracksNew(candidateTracksTemp, cells[itCell],
                      usedCells2);  // Move back to using used cells here? With low chi2/ndof?

      for (int itTrack = 0; itTrack < candidateTracksTemp.size(); itTrack++) {
        if (candidateTracksTemp[itTrack]->size() >= (m_minClustersOnTrack - 1))
          candidateTracks.push_back(candidateTracksTemp[itTrack]);
        else
          delete candidateTracksTemp[itTrack];
      }

      // Debug plotting
      if (m_debugPlots && m_eventNumber == 2) {
        m_canvConformalEventDisplayAcceptedCells->cd();
        for (int iTr = 0; iTr < candidateTracks.size(); iTr++) {
          cellularTrack* track = candidateTracks[iTr];
          for (unsigned int trackCell = 0; trackCell < track->size(); trackCell++) {
            drawline((*track)[trackCell]->getStart(), (*track)[trackCell]->getEnd(), track->size() - trackCell);
          }
        }
      }

      // Look at the candidate tracks and fit them (+pick the best chi2/ndof)
      //       streamlog_out(DEBUG)<<"- produced "<<candidateTracks.size()<<" candidate tracks"<<std::endl;
      if (debugSeed && kdhit == debugSeed)
        streamlog_out(DEBUG) << "- produced " << candidateTracks.size() << " candidate tracks" << std::endl;
      if (candidateTracks.size() == 0)
        continue;
      std::vector<double>   chi2ndof;
      std::vector<KDTrack*> bestTracks;

      getFittedTracks(bestTracks, candidateTracks,
                      usedCells);  // Returns all tracks at the moment, not lowest chi2 CHANGE ME

      // Store track(s) for later
      cellTracks.insert(cellTracks.end(), bestTracks.begin(), bestTracks.end());
    }

    // All tracks leading back to the seed hit have now been found. Decide which are the feasible candidates (may be more than 1)
    if (debugSeed && kdhit == debugSeed)
      streamlog_out(DEBUG) << "== final number of candidate tracks to this seed hit: " << cellTracks.size() << std::endl;
    if (cellTracks.size() == 0) {
      // Clean up
      for (unsigned int itCell = 0; itCell < cells.size(); itCell++)
        delete cells[itCell];
      continue;
    }

    std::vector<KDTrack*> bestTracks = cellTracks;  //CHANGE ME - temp to give all tracks
    //      std::vector<KDTrack*> bestTracks = getLowestChi2(cellTracks,cellTracksChi2ndof,chi2ndof);
    //     streamlog_out(DEBUG)<<"== final number of stored tracks to this seed hit: "<<bestTracks.size()<<std::endl;
    if (debugSeed && kdhit == debugSeed) {
      streamlog_out(DEBUG) << "== final number of stored tracks to this seed hit: " << bestTracks.size() << std::endl;
      for (int itBest = 0; itBest < bestTracks.size(); itBest++)
        streamlog_out(DEBUG) << "- track " << itBest << " has chi2/ndof " << bestTracks[itBest]->chi2ndof() << std::endl;
    }

    // Could now think to do the full helix fit and apply a chi2 cut. TODO

    // Store the final CA tracks. First sort them by length, so that if we produce a clone with just 1 hit missing, the larger track is taken
    std::sort(bestTracks.begin(), bestTracks.end(), sort_by_length);
    for (unsigned int itTrack = 0; itTrack < bestTracks.size(); itTrack++) {
      bool bestTrackUsed = false;

      // Cut on chi2
      //        if(collection != trackerHitCollections.size() && (bestTracks[itTrack]->chi2ndof() > m_chi2cut || bestTracks[itTrack]->chi2ndofZS() > m_chi2cut)) {
      double chi2cut = m_chi2cut;
      if (bestTracks[itTrack]->pt() < 5.)
        chi2cut = 1000.;

      if (bestTracks[itTrack]->chi2ndof() > chi2cut || bestTracks[itTrack]->chi2ndofZS() > chi2cut) {
        delete bestTracks[itTrack];
        bestTracks[itTrack] = NULL;
        continue;
      }

      // Check if the new track is a clone
      bool clone = false;

      for (int existingTrack = 0; existingTrack < conformalTracks.size(); existingTrack++) {
        const int nOverlappingHits = overlappingHits(bestTracks[itTrack], conformalTracks[existingTrack]);
        if (nOverlappingHits >= 2) {
          clone = true;

          // Calculate the new and existing chi2 values
          double newchi2 = (bestTracks[itTrack]->chi2ndofZS() * bestTracks[itTrack]->chi2ndofZS() +
                            bestTracks[itTrack]->chi2ndof() * bestTracks[itTrack]->chi2ndof());
          double oldchi2 = (conformalTracks[existingTrack]->chi2ndofZS() * conformalTracks[existingTrack]->chi2ndofZS() +
                            conformalTracks[existingTrack]->chi2ndof() * conformalTracks[existingTrack]->chi2ndof());

          double deltachi2ZS = (bestTracks[itTrack]->chi2ndofZS() - conformalTracks[existingTrack]->chi2ndofZS());
          double deltachi2   = (bestTracks[itTrack]->chi2ndof() - conformalTracks[existingTrack]->chi2ndof());

          // If the new track is an existing track + segment, make sure the increase in chi2 is not too much
          if (nOverlappingHits == conformalTracks[existingTrack]->m_clusters.size()) {
            // Increase in chi2 is too much (double)
            if ((newchi2 - oldchi2) > oldchi2)
              break;

            // Otherwise replace the existing track
            delete conformalTracks[existingTrack];
            conformalTracks[existingTrack] = bestTracks[itTrack];
            bestTrackUsed                  = true;
          }

          // If the new track is a subtrack of an existing track, don't consider it further (already try removing bad hits from tracks
          else if (nOverlappingHits == bestTracks[itTrack]->m_clusters.size())
            break;

          // Otherwise take the longest if the delta chi2 is not too much
          else if (bestTracks[itTrack]->m_clusters.size() >=
                   conformalTracks[existingTrack]->m_clusters.size()) {  // New track longer/equal in length

            // Increase in chi2 is too much (double)
            if ((newchi2 - oldchi2) > oldchi2)
              break;

            // Otherwise take it
            delete conformalTracks[existingTrack];
            conformalTracks[existingTrack] = bestTracks[itTrack];
            bestTrackUsed                  = true;
          } else if (bestTracks[itTrack]->m_clusters.size() <
                     conformalTracks[existingTrack]->m_clusters.size()) {  // Old track longer

            // Must improve chi2 by factor two
            if ((newchi2 - 0.5 * oldchi2) > 0.)
              break;

            // Otherwise take it
            delete conformalTracks[existingTrack];
            conformalTracks[existingTrack] = bestTracks[itTrack];
            bestTrackUsed                  = true;
          }

          break;
        }
      }

      // If not a clone, save the new track
      if (!clone) {
        conformalTracks.push_back(bestTracks[itTrack]);
        bestTrackUsed = true;
        if (debugSeed && kdhit == debugSeed)
          streamlog_out(DEBUG) << "== Pushing back best track with chi2/ndof " << bestTracks[itTrack]->chi2ndof()
                               << std::endl;
        else
          streamlog_out(DEBUG) << "Pushing back best track with chi2/ndof " << bestTracks[itTrack]->chi2ndof() << std::endl;
      }

      if (not bestTrackUsed) {
        delete bestTracks[itTrack];
        bestTracks[itTrack] = NULL;
      }

    }  // end for besttracks

    // Clean up
    for (unsigned int itCell = 0; itCell < cells.size(); itCell++)
      delete cells[itCell];
  }
}

// Take a collection of tracks and try to extend them into the collection of clusters passed.
void ConformalTracking::extendTracks(std::vector<KDTrack*>& conformalTracks, std::vector<KDCluster*>& collection,
                                     KDTree* nearestNeighbours) {
  // Loop over all current tracks. At the moment this is a "stupid" algorithm: it will simply try to add every
  // hit in the collection to every track, and keep the ones thta have a good chi2. In fact, it will extrapolate
  // the track and do a nearest neighbours search, but this seemed to fail for some reason, TODO!

  streamlog_out(DEBUG) << "EXTENDING tracks" << std::endl;
  if (collection.size() == 0)
    return;
  int nTracks = conformalTracks.size();
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

    // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
    Cell* seedCell = new Cell(conformalTracks[currentTrack]->m_clusters[1], conformalTracks[currentTrack]->m_clusters[0]);

    // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
    double     searchDistance = m_maxDistance;
    KDCluster* fakeHit = extrapolateCell(seedCell, searchDistance / 2.);  // TODO: make this search a function of radius
    VecCluster results2;
    nearestNeighbours->allNeighboursInRadius(fakeHit, 10 * searchDistance / 2., results2);
    std::sort(results2.begin(), results2.end(), sort_by_radiusKD);
    delete fakeHit;

    // Loop over all hits found and check if any have a sufficiently good delta Chi2
    vector<KDCluster*> goodHits;
    KDCluster*         bestCluster = NULL;
    double             bestChi2    = 0.;

    // Sort the hit collection by layer
    std::sort(collection.begin(), collection.end(), sort_by_layer);

    unsigned int nKDHits = collection.size();
    //streamlog_out(DEBUG)<<"STARTING"<<std::endl;
    //      for(int newHit=0;newHit<results2.size();newHit++){
    for (unsigned int nKDHit = 0; nKDHit < nKDHits; nKDHit++) {
      // Get the kdHit and check if it has already been used (assigned to a track)
      KDCluster* kdhit = collection[nKDHit];
      //streamlog_out(DEBUG)<<"Detector "<<kdhit->getSubdetector()<<", layer "<<kdhit->getLayer()<<", side "<<kdhit->getSide()<<std::endl;

      // If this hit is on a new layer, then add the hit from the previous layer and start anew
      if (bestCluster != NULL && !(kdhit->sameLayer(bestCluster))) {
        bestCluster->used(true);
        conformalTracks[currentTrack]->add(bestCluster);
        conformalTracks[currentTrack]->linearRegression();
        conformalTracks[currentTrack]->linearRegressionConformal();
        bestCluster = NULL;
      }

      if (kdhit->used())
        continue;

      // First check that the hit is not wildly away from the track (make cell and check angle)
      //        Cell* extensionCell = new Cell(conformalTracks[currentTrack]->m_clusters[0],results2[newHit]);
      Cell*  extensionCell = new Cell(conformalTracks[currentTrack]->m_clusters[0], kdhit);
      double cellAngle     = seedCell->getAngle(extensionCell);
      double cellAngleRZ   = seedCell->getAngleRZ(extensionCell);
      delete extensionCell;
      if (cellAngle > 3. * m_maxCellAngle || cellAngleRZ > 3. * m_maxCellAngleRZ) {
        continue;
      }

      // Now fit the track with the new hit and check the increase in chi2
      double deltaChi2 = fitWithPoint(*(conformalTracks[currentTrack]),
                                      kdhit);  //conformalTracks[currentTrack]->deltaChi2(results2[newHit]);

      // We have an estimate of the pT here, could use it in the chi2 criteria
      //          if(deltaChi2 > (10.*1000./conformalTracks[currentTrack]->m_pT)) continue;
      if (deltaChi2 > 1000.)
        continue;

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

    delete seedCell;

  }  // End of loop over tracks
}

// Extend seed cells
void ConformalTracking::extendSeedCells(std::vector<Cell*>& cells, KDTree* nearestNeighbours, bool extendingTrack,
                                        std::vector<KDCluster*> debugHits) {
  unsigned int nCells   = 0;
  int          depth    = 0;
  int          startPos = 0;

  // Keep track of existing cells in case there are branches in the track
  std::map<KDCluster*, std::vector<Cell*>> existingCells;

  // Try to create all "downstream" cells until no more can be added
  while (cells.size() != nCells) {
    // Extend all cells with depth N. In the next iteration, look at cells with depth N+1
    nCells = cells.size();
    for (unsigned int itCell = startPos; itCell < nCells; itCell++) {
      // Get the end point of the cell (to search for neighbouring hits to form new cells connected to this one)
      KDCluster* hit            = cells[itCell]->getEnd();
      double     searchDistance = m_maxDistance;  //hit->getR();
                                                  //      if(searchDistance > hit->getR()) searchDistance = 1.2*hit->getR();

      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      KDCluster* fakeHit =
          extrapolateCell(cells[itCell], searchDistance / 2.);  // TODO: make this search a function of radius
      VecCluster results;
      nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25 * searchDistance / 2., results);
      delete fakeHit;

      if (extendingTrack)
        streamlog_out(DEBUG) << "Found " << results.size() << " neighbours from cell extrapolation" << std::endl;
      // Make new cells pointing inwards
      for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
        if (extendingTrack)
          streamlog_out(DEBUG) << "looking at neighbour " << neighbour << std::endl;
        // Get the neighbouring hit
        KDCluster* nhit = results[neighbour];

        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
        //        if(used.count(nhit)){if(extendingTrack)streamlog_out(DEBUG)<<"- used"<<std::endl; continue;}
        if (nhit->used())
          continue;
        if (hit->sameLayer(nhit)) {
          if (extendingTrack)
            streamlog_out(DEBUG) << "- same layer" << std::endl;
          continue;
        }
        if (nhit->getR() >= hit->getR()) {
          if (extendingTrack)
            streamlog_out(DEBUG) << "- higher radius" << std::endl;
          continue;
        }

        // Check if this cell already exists (rejoining branch) FIXME - allows rejoining a branch without checking cell angles
        if (existingCells.count(hit) != 0) {
          bool alreadyExists  = false;
          int  nExistingCells = existingCells[hit].size();
          for (int iterExisting = 0; iterExisting < nExistingCells; iterExisting++) {
            if (existingCells[hit][iterExisting]->getEnd() == nhit) {
              alreadyExists = true;

              // Check if cell angle is too large to rejoin
              if (cells[itCell]->getAngle(existingCells[hit][iterExisting]) > m_maxCellAngle ||
                  cells[itCell]->getAngleRZ(existingCells[hit][iterExisting]) > m_maxCellAngleRZ)
                continue;

              // Otherwise add the path
              cells[itCell]->setTo(existingCells[hit][iterExisting]);
              existingCells[hit][iterExisting]->setFrom(cells[itCell]);
              updateCell(existingCells[hit][iterExisting]);
            }
          }
          if (alreadyExists)
            continue;
        }

        // Make the new cell
        Cell* cell = new Cell(hit, nhit);

        // Check if the new cell is compatible with the previous cell (angle between the two is acceptable)
        //        if( cells[itCell]->getAngle(cell) > (m_maxCellAngle*exp(-0.001/nhit->getR())) ){
        if (cells[itCell]->getAngle(cell) > m_maxCellAngle || cells[itCell]->getAngleRZ(cell) > m_maxCellAngleRZ) {
          // Debug plotting
          //          if(m_debugPlots && m_eventNumber == 0){
          //            m_canvConformalEventDisplayAllCells->cd();
          //            drawline(hit,nhit,cells[itCell]->getWeight()+2,3);
          //          }

          delete cell;
          continue;
        }

        // Set the information about which cell this new cell is attached to and store it
        cell->setFrom(cells[itCell]);
        cells[itCell]->setTo(cell);
        cells.push_back(cell);
        existingCells[hit].push_back(cell);

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

//===================================
// Cellular Track Functions
//===================================

// New test at creating cellular tracks. In this variant, don't worry about clones etc, give all possible routes back to the seed cell. Then cut
// on number of clusters on each track, and pass back (good tracks to then be decided based on best chi2
void ConformalTracking::createTracksNew(std::vector<cellularTrack*>& finalcellularTracks, Cell* seedCell,
                                        std::map<Cell*, bool>&                                  usedCells) {
  // Final container to be returned
  std::vector<cellularTrack*> cellularTracks;

  // Make the first cellular track using the seed cell
  cellularTrack* seedTrack = new cellularTrack();
  cellularTracks.push_back(seedTrack);
  seedTrack->push_back(seedCell);

  // Now start to follow all paths back from this seed cell
  // While there are still tracks that are not finished (last cell weight 0), keep following their path
  while (toBeUpdated(cellularTracks)) {
    //   streamlog_out(DEBUG)<<"== Updating "<<cellularTracks.size()<<" tracks"<<std::endl;
    // Loop over all (currently existing) tracks
    int nTracks = cellularTracks.size();
    for (int itTrack = 0; itTrack < nTracks; itTrack++) {
      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom()->size() == 0) {
        //       streamlog_out(DEBUG)<<"-- Track "<<itTrack<<" is finished"<<std::endl;
        continue;
      }

      // While there is only one path leading from this cell, follow that path
      Cell* cell = cellularTracks[itTrack]->back();
      //     streamlog_out(DEBUG)<<"-- Track "<<itTrack<<" has "<<(*(cell->getFrom())).size()<<" cells attached to the end of it"<<std::endl;
      //      while(cell->getWeight() > 0 && (*(cell->getFrom())).size() == 1){
      while ((*(cell->getFrom())).size() == 1) {
        //       streamlog_out(DEBUG)<<"- simple extension"<<std::endl;
        // Get the cell that it attaches to
        Cell* parentCell = (*(cell->getFrom()))[0];
        // Attach it to the track and continue
        cellularTracks[itTrack]->push_back(parentCell);
        cell = parentCell;
      }

      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom()->size() == 0)
        continue;

      // If the weight is != 0 and there is more than one path to follow, branch the track (create a new one for each path)
      int nBranches = (*(cell->getFrom())).size();

      //     streamlog_out(DEBUG)<<"- making "<<nBranches<<" branches"<<std::endl;

      // For each additional branch make a new track
      for (int itBranch = 1; itBranch < nBranches; itBranch++) {
        cellularTrack* branchedTrack = new cellularTrack();
        (*branchedTrack)             = (*cellularTracks[itTrack]);
        cellularTracks.push_back(branchedTrack);
        branchedTrack->push_back((*(cell->getFrom()))[itBranch]);
      }

      // Keep the existing track for the first branch
      cellularTracks[itTrack]->push_back((*(cell->getFrom()))[0]);
    }
  }

  int nTracks = cellularTracks.size();
  for (int itTrack = 0; itTrack < nTracks; itTrack++) {
    //    if(cellularTracks[itTrack]->size() >= (m_minClustersOnTrack-1) )
    finalcellularTracks.push_back(cellularTracks[itTrack]);
  }

  return;
}

// Check if any of the tracks in a collection still have to be updated
bool ConformalTracking::toBeUpdated(std::vector<cellularTrack*> const& cellularTracks) {
  for (unsigned int iTrack = 0; iTrack < cellularTracks.size(); iTrack++)
    if (cellularTracks[iTrack]->back()->getFrom()->size() > 0) {
      return true;
    }
  return false;
}

// Given a list of connected cells (so-called cellular tracks), return the candidate(s) with lowest chi2/degrees of freedom.
// Attempt to remove hits to see if there is a significant improvement in the chi2/ndof, to protect against noise hits causing
// a good track to be discarded. If several candidates have low chi2/ndof and are not clones (limited sharing of hits) then
// return all of them. Given that there is no material scattering taken into account this helps retain low pt tracks, which
// may have worse chi2/ndof than ghosts/real tracks with an additional unrelated hit from the low pt track.
void ConformalTracking::getFittedTracks(std::vector<KDTrack*>& finalTracks, std::vector<cellularTrack*>& candidateTracks,
                                        std::map<Cell*, bool>&                                           usedCells) {
  // Make a container for all tracks being considered, initialise variables
  std::vector<KDTrack*> trackContainer;
  //  std::vector<double> trackChi2ndofs;

  // Loop over all candidate tracks and do an inital fit to get the track angle (needed to calculate the
  // hit errors for the error-weighted fit)
  for (unsigned int itTrack = 0; itTrack < candidateTracks.size(); itTrack++) {
    // If there are not enough hits on the track, ignore it
    if (candidateTracks[itTrack]->size() < (m_minClustersOnTrack - 2)) {
      delete candidateTracks[itTrack];
      candidateTracks[itTrack] = NULL;
      continue;
    }

    // Make the fitting object. TGraphErrors used for 2D error-weighted fitting
    KDTrack* track = new KDTrack();

    // Loop over all hits and add them to the fitter (and track)
    double     npoints = 0.;
    KDCluster* kdStart = (*candidateTracks[itTrack])[0]->getEnd();
    track->add(kdStart);
    npoints++;

    for (unsigned int trackCell = 0; trackCell < (*candidateTracks[itTrack]).size(); trackCell++) {
      KDCluster* kdEnd = (*candidateTracks[itTrack])[trackCell]->getStart();
      track->add(kdEnd);
      npoints++;
    }
    track->linearRegression();
    track->linearRegressionConformal();               //FCC study
    double chi2ndof = track->chi2() / (npoints - 2);  //FCC study

    // We try to see if there are spurious hits causing the chi2 to be very large. This would cause us to throw away
    // good tracks with perhaps just a single bad hit. Try to remove each hit and see if fitting without it causes a
    // significant improvement in chi2/ndof

    // Loop over each hit (starting at the back, since we will use the 'erase' function to get rid of them)
    // and see if removing it improves the chi2/ndof
    double removed = 0;
    if (chi2ndof > m_chi2cut &&
        chi2ndof <
            m_chi2cut) {  //CHANGE ME?? Upper limit to get rid of really terrible tracks (temp lower changed from 0 to m_chi2cut)
      for (int point = npoints - 1; point >= 0; point--) {
        // Stop if we would remove too many points on the track to meet the minimum hit requirement (or if the track has more
        // than 2 hits removed)
        if ((npoints - removed - 1) < m_minClustersOnTrack || removed == 2)
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
    trackContainer.push_back(track);
    //    trackChi2ndofs.push_back(chi2ndof);

    delete candidateTracks[itTrack];

  }  // end for candidateTracks

  // Now have all sets of conformal tracks and their chi2/ndof. Decide which tracks to send back, ie. the one with
  // lowest chi2/ndof, and possibly others if they are not clones and have similar chi2 value
  getLowestChi2(finalTracks, trackContainer);

  // Send back the final set of tracks
  return;
}

// Pick the lowest chi2/ndof KDTrack from a list of possible tracks, and additionally return other tracks in the collection with similar chi2/ndof values that don't share many hits
void ConformalTracking::getLowestChi2(std::vector<KDTrack*>& finalTracks, std::vector<KDTrack*> trackContainer) {
  // Get the lowest chi2/ndof value from the given tracks
  //  double lowestChi2ndof = *std::min_element(trackChi2ndofs.begin(),trackChi2ndofs.end());
  KDTrack* lowestChi2ndofTrack = trackContainer[0];
  double   lowestChi2ndof      = lowestChi2ndofTrack->chi2ndof();

  for (unsigned int itTrack = 0; itTrack < trackContainer.size(); itTrack++) {
    if (trackContainer[itTrack]->chi2ndof() < lowestChi2ndof) {
      lowestChi2ndof      = trackContainer[itTrack]->chi2ndof();
      lowestChi2ndofTrack = trackContainer[itTrack];
    }
  }

  //finalTracks.push_back(lowestChi2ndofTrack);
  //return;

  // Loop over all other tracks and decide whether or not to save them
  for (unsigned int itTrack = 0; itTrack < trackContainer.size(); itTrack++) {
    // Look at the difference in chi2/ndof - we want to keep tracks with similar chi2/ndof. If they
    // are clones then take the longest
    if ((trackContainer[itTrack]->chi2ndof() - lowestChi2ndof) < 10.) {
      // Store this track
      finalTracks.push_back(trackContainer[itTrack]);

    } else {
      delete trackContainer[itTrack];
      trackContainer[itTrack] = NULL;
    }
  }

  return;
}

void ConformalTracking::updateCell(Cell* cell) {
  if ((*(cell->getTo())).size() != 0) {
    for (unsigned int i = 0; i < (*(cell->getTo())).size(); i++) {
      (*(cell->getTo()))[i]->update(cell);
      updateCell((*(cell->getTo()))[i]);
    }
  }
}

// Function to extrapolate along a cell in conformal space, producing a fake hit
// a given distance away from the cell endpoint
KDCluster* ConformalTracking::extrapolateCell(Cell* cell, double distance) {
  // Fake cluster to be returned
  KDCluster* extrapolatedCluster = new KDCluster();

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
void ConformalTracking::extendTrack(KDTrack* track, std::vector<cellularTrack*> trackSegments,
                                    std::map<KDCluster*, bool>& used, std::map<Cell*, bool>& usedCells) {
  //  cout<<"== extending track, have "<<trackSegments.size()<<" candidates"<<endl;
  // For each track segment, perform a kalman filter on the addition points and chose the track extension with the
  // best delta chi2.
  double             bestChi2 = 1.e9;
  double             bestNpoints;
  vector<KDCluster*> finalHits;
  double             bestChi2Conformal, bestChi2ZS;

  //  for(int i=0;i<track->clusters().size();i++)streamlog_out(DEBUG)<<"- cluster "<<i<<" has u,v "<<track->clusters()[i]->getU()<<","<<track->clusters()[i]->getV()<<std::endl;

  for (int nTrackExtension = 0; nTrackExtension < trackSegments.size(); nTrackExtension++) {
    //   streamlog_out(DEBUG)<<"Track extension "<<nTrackExtension<<std::endl;
    vector<KDCluster*> newHits;

    // Add each of the new clusters
    KDCluster* kdStart = (*trackSegments[nTrackExtension])[0]->getEnd();
    newHits.push_back(kdStart);
    //   streamlog_out(DEBUG)<<"- cluster "<<" has u,v "<<kdStart->getU()<<","<<kdStart->getV()<<std::endl;

    for (unsigned int trackCell = 0; trackCell < trackSegments[nTrackExtension]->size() - 2; trackCell++) {
      KDCluster* kdEnd = (*trackSegments[nTrackExtension])[trackCell]->getStart();
      newHits.push_back(kdEnd);
      //     streamlog_out(DEBUG)<<"- cluster "<<" has u,v "<<kdEnd->getU()<<","<<kdEnd->getV()<<std::endl;
    }

    double newChi2Conformal, newChi2ZS;
    double newChi2 = fitWithExtension(*track, newHits, newChi2Conformal, newChi2ZS);

    if (newChi2 < bestChi2) {
      finalHits.clear();
      for (int i = 0; i < newHits.size(); i++)
        finalHits.push_back(newHits[i]);
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
    //   streamlog_out(DEBUG)<<"- new chi2/ndof: "<<track->chi2ndof()<<", chi2ZS/ndof: "<<track->chi2ndofZS()<<std::endl;
  }
}

//===================================
// Fitting Functions
//===================================

double ConformalTracking::fitWithExtension(KDTrack track, std::vector<KDCluster*> hits, double& newChi2, double& newChi2ZS) {
  // Add the point to the track
  for (int i = 0; i < hits.size(); i++)
    track.add(hits[i]);
  int npoints = track.nPoints();

  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal();

  double chi2ndof = sqrt(track.chi2ndofZS() * track.chi2ndofZS() + track.chi2ndof() * track.chi2ndof());
  newChi2         = track.chi2ndof();
  newChi2ZS       = track.chi2ndofZS();
  return chi2ndof;
}

// Add a point to a track and return the delta chi2
double ConformalTracking::fitWithPoint(KDTrack kdTrack, KDCluster* hit) {
  double chi2   = kdTrack.chi2ndof();
  double chi2zs = kdTrack.chi2ndofZS();
  kdTrack.add(hit);
  kdTrack.linearRegression();
  kdTrack.linearRegressionConformal();
  double newchi2   = kdTrack.chi2ndof();
  double newchi2zs = kdTrack.chi2ndofZS();

  return (newchi2 - chi2 + newchi2zs - chi2zs);
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
void ConformalTracking::getCollection(LCCollection*& collection, std::string collectionName, LCEvent* evt) {
  try {
    collection = evt->getCollection(collectionName);
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

// Function to check if two KDtracks contain several hits in common
int ConformalTracking::overlappingHits(const KDTrack* track1, const KDTrack* track2) {
  int nHitsInCommon = 0;
  for (int hit = 0; hit < track1->m_clusters.size(); hit++) {
    if (std::find(track2->m_clusters.begin(), track2->m_clusters.end(), track1->m_clusters[hit]) != track2->m_clusters.end())
      nHitsInCommon++;
  }
  return nHitsInCommon;
}

//===================================
// Debug Functions for Reconstruction
//===================================

// Debug function - checks if a track will be associated to an MC particle or not
double ConformalTracking::checkReal(KDTrack* track, std::map<KDCluster*, MCParticle*> kdParticles,
                                    std::map<MCParticle*, bool>&                   reconstructed,
                                    std::map<MCParticle*, std::vector<KDCluster*>> MCparticleHits) {
  // Store all mcparticles associated to this track
  std::vector<MCParticle*> particles;
  std::map<MCParticle*, double> particleHits;
  double nHits = 0.;

  // Get the clusters from this track
  std::vector<KDCluster*> clusters = track->m_clusters;

  // Loop over all hits and see which particle they are associated to
  for (int itCluster = 0; itCluster < clusters.size(); itCluster++) {
    // Get the hit
    KDCluster* cluster = clusters[itCluster];
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
  for (int iPart = 0; iPart < particles.size(); iPart++) {
    if (particleHits[particles[iPart]] > bestHits) {
      bestHits     = particleHits[particles[iPart]];
      bestParticle = particles[iPart];
    }
  }

  // Calculate the purity
  double purity = bestHits / nHits;
  streamlog_out(DEBUG) << "Number of hits on track: " << nHits << ". Good hits: " << bestHits << ". Purity: " << purity
                       << ". Pt: " << sqrt(bestParticle->getMomentum()[0] * bestParticle->getMomentum()[0] +
                                           bestParticle->getMomentum()[1] * bestParticle->getMomentum()[1])
                       << ". Track chi2/ndof: " << track->chi2ndof() << ". Chi2/ndof in SZ fit: " << track->chi2ndofZS()
                       << std::endl;

  // Check the track pt estimate
  TLorentzVector mc_helper;
  mc_helper.SetPxPyPzE(bestParticle->getMomentum()[0], bestParticle->getMomentum()[1], bestParticle->getMomentum()[2],
                       bestParticle->getEnergy());

  std::cout << "Particle at theta = " << mc_helper.Theta() * 180. / M_PI << ", with pt = " << mc_helper.Pt()
            << ". pt estimate: " << track->pt() << std::endl;

  // Check if any hits are missing
  std::vector<KDCluster*> mcHits     = MCparticleHits[bestParticle];
  int                     uniqueHits = getUniqueHits(mcHits);

  streamlog_out(DEBUG) << "Track should contain " << uniqueHits << " hits, "
                       << ((uniqueHits > bestHits) ? "is missing hits" : "all hits found.") << std::endl;

  std::vector<KDCluster*> trackHits = track->m_clusters;
  for (int th = 0; th < trackHits.size(); th++) {
    streamlog_out(DEBUG) << "Hit " << th << " u = " << trackHits[th]->getU() << " v = " << trackHits[th]->getV()
                         << " x = " << trackHits[th]->getX() << " y = " << trackHits[th]->getY()
                         << " z = " << trackHits[th]->getZ() << std::endl;

    // Check which particle the hit belongs to
    MCParticle* particle  = kdParticles[trackHits[th]];
    double      mcVertexX = particle->getVertex()[0];
    double      mcVertexY = particle->getVertex()[1];
    double      mcVertexR = sqrt(pow(mcVertexX, 2) + pow(mcVertexY, 2));

    streamlog_out(DEBUG) << "Come from particle with id " << particle->getPDG() << " produced at vertexR = " << mcVertexR;

    if (particle->getParents().size() != 0) {
      streamlog_out(DEBUG) << ". Comes from particle with id " << particle->getParents()[0]->getPDG();
    }
    streamlog_out(DEBUG) << std::endl;
  }

  //  for(int itCluster=0;itCluster<clusters.size();itCluster++)streamlog_out(DEBUG)<<"Hit "<<itCluster<<" has position: "<<clusters[itCluster]->getU()<<","<<clusters[itCluster]->getV()<<std::endl;
  streamlog_out(DEBUG) << "== Terms in conformal fit: " << track->m_intercept << ", " << track->m_gradient << ", "
                       << track->m_quadratic << std::endl;
  streamlog_out(DEBUG) << "== Terms in zs fit: " << track->m_interceptZS << ", " << track->m_gradientZS << std::endl;

  track->linearRegressionConformal(true);
  track->calculateChi2SZ(NULL, true);
  if (purity >= m_purity) {
    reconstructed[bestParticle] = true;
  }

  // Return the purity
  return purity;
}

// Debug function - gets number of unique hits
int ConformalTracking::getUniqueHits(std::vector<KDCluster*> hits) {
  int nUniqueHits = 0;

  for (int iHit = 0; iHit < hits.size(); iHit++) {
    // Get the cluster
    KDCluster* hit    = hits[iHit];
    bool       sameID = false;

    // Check if any following clusters have the same detector id
    for (int iHitLater = iHit + 1; iHitLater < hits.size(); iHitLater++) {
      // Get the new cluster
      KDCluster* nhit = hits[iHitLater];

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
                                                   std::map<MCParticle*, std::vector<KDCluster*>> particleHits,
                                                   KDTree* nearestNeighbours) {
  // Get the hits for this MC particle
  std::vector<KDCluster*> trackHits = particleHits[particle];

  // Sort the hits from larger to smaller radius
  std::sort(trackHits.begin(), trackHits.end(), sort_by_radiusKD);
  streamlog_out(DEBUG) << "Track contains " << getUniqueHits(trackHits) << " unique hits. Track size is " << trackHits.size()
                       << std::endl;

  for (int th = 0; th < trackHits.size(); th++)
    streamlog_out(DEBUG) << "Hit " << th << " u = " << trackHits[th]->getU() << " v = " << trackHits[th]->getV()
                         << " x = " << trackHits[th]->getX() << " y = " << trackHits[th]->getY()
                         << " z = " << trackHits[th]->getZ() << ". " << (trackHits[th]->used() ? "Used!" : "") << std::endl;

  // Take the seed hit and build cell to 2nd
  KDCluster* seedHit = trackHits[0];
  VecCluster results;
  streamlog_out(DEBUG) << "- getting nearest neighbours for seed" << std::endl;
  //  nearestNeighbours->allNeighboursInRadius(seedHit, m_maxDistance, results);
  double theta = seedHit->getTheta();
  nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);
  streamlog_out(DEBUG) << "- got nearest neighbours for seed" << std::endl;

  // If we don't produce the seed cell we have failed
  if (std::find(results.begin(), results.end(), trackHits[1]) == results.end()) {
    streamlog_out(DEBUG) << "- seed cell not produced, neighbour not inside search window" << std::endl;
    double deltaU = fabs(trackHits[1]->getU() - trackHits[0]->getU());
    double deltaV = fabs(trackHits[1]->getV() - trackHits[0]->getV());
    streamlog_out(DEBUG) << "- distance between hits is " << sqrt(deltaU * deltaU + deltaV * deltaV)
                         << " and search distance is " << m_maxDistance << std::endl;
    streamlog_out(DEBUG) << "- theta of hits is " << trackHits[0]->getTheta() << " and " << trackHits[1]->getTheta()
                         << ", delta theta = " << trackHits[0]->getTheta() - trackHits[1]->getTheta() << std::endl;
    //    return;
  }

  if (trackHits[1]->getR() >= trackHits[0]->getR()) {
    streamlog_out(DEBUG) << "- seed cell not produced, neighbour at higher radius" << std::endl;
    //    return;
  }

  double length = sqrt((trackHits[1]->getU() - trackHits[0]->getU()) * (trackHits[1]->getU() - trackHits[0]->getU()) +
                       (trackHits[1]->getV() - trackHits[0]->getV()) * (trackHits[1]->getV() - trackHits[0]->getV()));
  if (length > m_maxDistance) {
    streamlog_out(DEBUG) << "- seed cell not produced, neighbour too far away" << std::endl;
    //    return;
  }

  // Make the seed cell and extrapolate until all hits found
  Cell*              seedCell = new Cell(trackHits[0], trackHits[1]);
  std::vector<Cell*> cells;
  cells.push_back(seedCell);
  streamlog_out(DEBUG) << "have seed cell" << std::endl;
  // CHECK CELL LENGTH!! TO DO

  int  hitNumber     = 1;  // Current hit after the seed cell has been formed
  bool trackBuilding = true;
  while (hitNumber < (trackHits.size() - 1)) {
    // Extend the current cell
    double searchDistance = m_maxDistance;  //hit->getR();
    //    if(searchDistance > trackHits[hitNumber]->getR()) searchDistance = trackHits[hitNumber]->getR();

    // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
    KDCluster* fakeHit = extrapolateCell(cells.back(), searchDistance / 2.);  // TODO: make this search a function of radius
    VecCluster results2;
    nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25 * searchDistance / 2., results2);
    delete fakeHit;

    // Check if we found the hit we wanted
    if (std::find(results2.begin(), results2.end(), trackHits[hitNumber + 1]) == results2.end()) {
      streamlog_out(DEBUG) << "- could not find hit number " << hitNumber + 1 << " inside the search window" << std::endl;
      double deltaU = fabs(trackHits[hitNumber + 1]->getU() - trackHits[hitNumber]->getU());
      double deltaV = fabs(trackHits[hitNumber + 1]->getV() - trackHits[hitNumber]->getV());
      streamlog_out(DEBUG) << "- distance between hits is " << sqrt(deltaU * deltaU + deltaV * deltaV)
                           << " and search distance is " << m_maxDistance << std::endl;
      //      for(int c=0;c<cells.size();c++) delete cells[c];
      //      return;
    }

    // Check radial conditions
    if (trackHits[hitNumber + 1]->getR() >= trackHits[hitNumber]->getR()) {
      streamlog_out(DEBUG) << "- cell " << hitNumber << " not produced, neighbour at higher radius" << std::endl;
      //      return;
    }

    // Make the cell for extrapolation in the next round
    Cell* cell = new Cell(trackHits[hitNumber], trackHits[hitNumber + 1]);
    streamlog_out(DEBUG) << "have new cell" << std::endl;

    // Check if the cell would be killed by cuts
    if (cells.back()->getAngle(cell) > m_maxCellAngle || cells.back()->getAngleRZ(cell) > m_maxCellAngleRZ) {
      if (cells.back()->getAngle(cell) > m_maxCellAngle)
        streamlog_out(DEBUG) << "- cell " << hitNumber << " killed by angular cut. Angle is " << cells.back()->getAngle(cell)
                             << std::endl;
      if (cells.back()->getAngleRZ(cell) > m_maxCellAngleRZ)
        streamlog_out(DEBUG) << "- cell " << hitNumber << " killed by RZ angular cut. Angle is "
                             << cells.back()->getAngleRZ(cell) << std::endl;
      //      for(int c=0;c<cells.size();c++) delete cells[c];
      //      return;
    }

    // Set the information about which cell this new cell is attached to and store it
    cell->setFrom(cells.back());
    cells.back()->setTo(cell);
    //    cell->setWeight(hitNumber);
    cells.push_back(cell);

    // Check if we have finished
    //    if(hitNumber == (trackHits.size()-1)) trackBuilding = false;
    hitNumber++;
  }

  // Now test the built-in functions to make CellularTracks from cells
  std::vector<cellularTrack*> cellularTracks;
  std::map<Cell*, bool> usedCells;
  createTracksNew(cellularTracks, cells.back(), usedCells);

  if (cellularTracks.size() != 1) {
    streamlog_out(DEBUG) << "- cellular track not produced from cells. Returned " << cellularTracks.size() << " candidates"
                         << std::endl;
    for (int c = 0; c < cells.size(); c++)
      delete cells[c];
    for (int ct = 0; ct < cellularTracks.size(); ct++)
      delete cellularTracks[ct];
    return;
  }

  // Check that all of the hits are on the track
  if (cellularTracks[0]->size() != (trackHits.size() - 1)) {
    streamlog_out(DEBUG) << "- failed to put all cells on track. Cellular track size: " << cellularTracks[0]->size()
                         << std::endl;
    streamlog_out(DEBUG) << "- number of cells held: " << cells.size() << std::endl;
    return;
  }

  // Now test the built-in functions to make KDTracks from CellularTracks
  std::vector<KDTrack*> finalTracks;
  getFittedTracks(finalTracks, cellularTracks, usedCells);

  if (finalTracks.size() != 1) {
    streamlog_out(DEBUG) << "- kd track not produced from cellular track. Returned " << finalTracks.size() << " candidates"
                         << std::endl;
    for (int c = 0; c < cells.size(); c++)
      delete cells[c];
    for (int t = 0; t < finalTracks.size(); t++)
      delete finalTracks[t];
    return;
  }

  // Test the chi2 criteria
  KDTrack* mcTrack = finalTracks[0];
  mcTrack->linearRegressionConformal(true);

  streamlog_out(DEBUG) << "== Track chi2/ndof is " << mcTrack->chi2ndof() << ", ZS chi2/ndof is " << mcTrack->chi2ndofZS()
                       << std::endl;
  streamlog_out(DEBUG) << "== Terms in conformal fit: " << mcTrack->m_intercept << ", " << mcTrack->m_gradient << ", "
                       << mcTrack->m_quadratic << std::endl;
  streamlog_out(DEBUG) << "== Terms in zs fit: " << mcTrack->m_interceptZS << ", " << mcTrack->m_gradientZS << std::endl;
  mcTrack->calculateChi2SZ(NULL, true);

  if (mcTrack->chi2ndof() > m_chi2cut || mcTrack->chi2ndofZS() > m_chi2cut) {
    if (mcTrack->chi2ndof() > m_chi2cut)
      streamlog_out(DEBUG) << "- track killed by chi2/ndof cut. Track chi2/ndof is " << mcTrack->chi2ndof() << std::endl;
    if (mcTrack->chi2ndofZS() > m_chi2cut)
      streamlog_out(DEBUG) << "- track killed by ZS chi2/ndof cut. Track chi2/ndof in ZS is " << mcTrack->chi2ndofZS()
                           << std::endl;
    for (int c = 0; c < cells.size(); c++)
      delete cells[c];
    for (int t = 0; t < finalTracks.size(); t++)
      delete finalTracks[t];
    return;
  }

  streamlog_out(DEBUG) << "== Should have produced this track" << std::endl;
  for (int c = 0; c < cells.size(); c++)
    delete cells[c];
  for (int t = 0; t < finalTracks.size(); t++)
    delete finalTracks[t];
  return;
}

// Draw a line on the current canvas
void ConformalTracking::drawline(KDCluster* hitStart, KDCluster* hitEnd, int colour, int style) {
  TLine* line = new TLine(hitStart->getU(), hitStart->getV(), hitEnd->getU(), hitEnd->getV());
  line->SetLineColor(colour);
  line->SetLineStyle(style);
  line->Draw();
}
