/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
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

#include "TCanvas.h"
#include "TFile.h"
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

  // Output collections - tracks
  registerOutputCollection(LCIO::TRACK, "SiTrackCollectionName", "Silicon track Collection Name", m_outputTrackCollection,
                           std::string("CATracks"));

  // Parameters for tracking
  registerProcessorParameter("DebugPlots", "Plots for debugging the tracking", m_debugPlots, bool(false));
  registerProcessorParameter("ThetaRange", "Angular range for initial cell seeding", m_thetaRange, double(0.4));
  registerProcessorParameter("MaxCellAngle", "Cut on angle between two cells for cell to be valid", m_maxCellAngle,
                             double(0.03));
  registerProcessorParameter("MaxDistance", "Maximum length of a cell (max. distance between two hits)", m_maxDistance,
                             double(0.02));
  registerProcessorParameter("MinClustersOnTrack", "Minimum number of clusters to create a track", m_minClustersOnTrack,
                             int(6));
}

bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2) {
  double radius1 =
      sqrt((hit1->getPosition()[0]) * (hit1->getPosition()[0]) + (hit1->getPosition()[1]) * (hit1->getPosition()[1]));
  double radius2 =
      sqrt((hit2->getPosition()[0]) * (hit2->getPosition()[0]) + (hit2->getPosition()[1]) * (hit2->getPosition()[1]));
  return (radius1 < radius2);
}

bool sort_by_radiusKD(KDCluster* hit1, KDCluster* hit2) {
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 > radius2);
}

bool sort_by_cellWeight(Cell* cell1, Cell* cell2) {
  int weight1 = cell1->getWeight();
  int weight2 = cell2->getWeight();
  return (weight1 > weight2);
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

  // Initialise histograms (if debug plotting on)
  if (m_debugPlots) {
    // Automatically add histograms to output root file
    AIDAProcessor::histogramFactory(this);

    // Histogram initailisation

    // Histograms for tuning parameters (cell angle cut, cell length cut)
    m_cellAngle        = new TH1F("cellAngle", "cellAngle", 1250, 0, 0.05);
    m_cellAngleRadius  = new TH2F("cellAngleRadius", "cellAngleRadius", 400, 0, 0.04, 1000, 0, 0.04);
    m_cellLengthRadius = new TH2F("cellLengthRadius", "cellLengthRadius", 300, 0, 0.03, 1000, 0, 0.04);
    m_cellAngleLength  = new TH2F("cellAngleLength", "cellAngleLength", 400, 0, 0.04, 300, 0, 0.03);

    m_cellAngleMC        = new TH1F("cellAngleMC", "cellAngleMC", 1250, 0, 0.05);
    m_cellAngleRadiusMC  = new TH2F("cellAngleRadiusMC", "cellAngleRadiusMC", 400, 0, 0.04, 1000, 0, 0.04);
    m_cellLengthRadiusMC = new TH2F("cellLengthRadiusMC", "cellLengthRadiusMC", 300, 0, 0.03, 1000, 0, 0.04);
    m_cellAngleLengthMC  = new TH2F("cellAngleLengthMC", "cellAngleLengthMC", 400, 0, 0.04, 300, 0, 0.03);

    // Histograms for "event display"
    m_conformalEvents       = new TH2F("conformalEvents", "conformalEvents", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    m_nonconformalEvents    = new TH2F("nonconformalEvents", "nonconformalEvents", 500, -1500, 1500, 500, -1500, 1500);
    m_conformalEventsRTheta = new TH2F("conformalEventsRTheta", "conformalEventsRTheta", 200, 0, 0.05, 632, -3.16, 3.16);
    m_conformalEventsMC     = new TH2F("conformalEventsMC", "conformalEventsMC", 1000, -0.05, 0.05, 1000, -0.05, 0.05);

    m_canvConformalEventDisplay   = new TCanvas("canvConformalEventDisplay", "canvConformalEventDisplay");
    m_canvConformalEventDisplayMC = new TCanvas("canvConformalEventDisplayMC", "canvConformalEventDisplayMC");
  }

  // Register this process
  Global::EVENTSEEDER->registerProcessor(this);
}

void ConformalTracking::processRunHeader(LCRunHeader* run) { ++m_runNumber; }

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

  // Object to store all of the track hit collections passed to the pattern recognition
  std::vector<LCCollection*> trackerHitCollections;

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
  }

  // Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec(LCIO::TRACK);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  /*
   Debug plotting. This section picks up tracks reconstructed using the cheated pattern recognition (TruthTrackFinder) and uses it to show
   values which are cut on during the tracking. This can be used to tune the cut ranges.
   */

  if (m_debugPlots) {
    // Draw the empty event display onto the canvas, so that cells can be added sequentially
    if (m_eventNumber == 0) {
      m_canvConformalEventDisplayMC->cd();
      m_conformalEventsMC->DrawCopy("");
    }

    // Get the collection of cheated tracks
    LCCollection* trueTrackCollection = 0;
    getCollection(trueTrackCollection, "SiTracks", evt);
    if (trackCollection == 0)
      return;

    // Loop over all tracks
    int nTracks = trueTrackCollection->getNumberOfElements();
    for (int itTrack = 0; itTrack < nTracks; itTrack++) {
      // Get the track
      Track* track = dynamic_cast<Track*>(trueTrackCollection->getElementAt(itTrack));

      // Get the hits
      const TrackerHitVec& hitVector = track->getTrackerHits();

      // Loop over hits and check the MC particle that they correspond to, and which subdetector they are on
      int nHits = hitVector.size();
      streamlog_out(DEBUG4) << "- New MC track with size " << nHits << std::endl;

      for (int itHit = 0; itHit < (nHits - 2); itHit++) {
        // Get the 3 tracker hits
        TrackerHitPlane* hit         = dynamic_cast<TrackerHitPlane*>(hitVector.at(itHit));
        TrackerHitPlane* nexthit     = dynamic_cast<TrackerHitPlane*>(hitVector.at(itHit + 1));
        TrackerHitPlane* nextnexthit = dynamic_cast<TrackerHitPlane*>(hitVector.at(itHit + 2));

        // Form the conformal clusters
        KDCluster* cluster0 = new KDCluster(hit);
        KDCluster* cluster1 = new KDCluster(nexthit);
        KDCluster* cluster2 = new KDCluster(nextnexthit);

        // Make the two cells connecting these three hits
        Cell* cell = new Cell(cluster0, cluster1);
        cell->setWeight(itHit);
        Cell* cell1 = new Cell(cluster1, cluster2);
        cell1->setWeight(itHit + 1);

        // Fill the debug/tuning plots
        double angleBetweenCells = cell->getAngle(cell1);
        double cell0Length = sqrt(pow(cluster0->getU() - cluster1->getU(), 2) + pow(cluster0->getV() - cluster1->getV(), 2));
        double cell1Length = sqrt(pow(cluster1->getU() - cluster2->getU(), 2) + pow(cluster1->getV() - cluster2->getV(), 2));

        m_cellAngleMC->Fill(angleBetweenCells);
        m_cellAngleRadiusMC->Fill(cluster2->getR(), angleBetweenCells);
        m_cellLengthRadiusMC->Fill(cluster0->getR(), cell0Length);
        m_cellAngleLengthMC->Fill(cell1Length, angleBetweenCells);

        // Draw cells on the first event
        if (m_eventNumber == 0) {
          // Fill the event display
          m_conformalEventsMC->Fill(cluster0->getU(), cluster0->getV());
          m_conformalEventsMC->Fill(cluster1->getU(), cluster1->getV());
          m_conformalEventsMC->Fill(cluster2->getU(), cluster2->getV());
          // Draw the cell lines on the event display
          m_canvConformalEventDisplayMC->cd();
          drawline(cluster0, cluster1, itHit + 1);
          drawline(cluster1, cluster2, itHit + 2);
        }
      }
    }

    // Draw the final set of conformal hits (on top of the cell lines)
    if (m_eventNumber == 0) {
      m_canvConformalEventDisplayMC->cd();
      m_conformalEventsMC->DrawCopy("same");
      // Draw the non-MC event display
      m_canvConformalEventDisplay->cd();
      m_conformalEvents->DrawCopy("");
    }
  }

  /*
   Now start doing things!
   */

  // Set up ID decoder
  UTIL::BitField64 m_encoder(lcio::ILDCellID0::encoder_string);

  // Some global containers to be used throughout the tracking. A collection of conformal hits will be made, with a link
  // pointing back to the corresponding cluster. A record of all used hits will be kept.
  std::map<KDCluster*, TrackerHitPlane*> kdClusterMap;
  std::vector<KDCluster*> kdClusters;
  std::map<KDCluster*, bool> used;
  std::vector<cellularTrack> finalCAtracks;

  // Loop over all input collections. Tracking will be attempted on the collection, then hits from the next collection
  // will be added to the unused hits already there.
  for (unsigned int collection = 0; collection < trackerHitCollections.size(); collection++) {
    // Loop over tracker hits and make conformal hit collection
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    for (int itHit = 0; itHit < nHits; itHit++) {
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(trackerHitCollections[collection]->getElementAt(itHit));

      // Make a new kd cluster
      KDCluster* kdhit = new KDCluster(hit);

      // Get subdetector information and add it to the kdhit
      const int celId = hit->getCellID0();
      m_encoder.setValue(celId);
      int subdet = m_encoder[lcio::ILDCellID0::subdet];
      int side   = m_encoder[lcio::ILDCellID0::side];
      int layer  = m_encoder[lcio::ILDCellID0::layer];
      kdhit->setDetectorInfo(subdet, side, layer);

      // Store the link between the two
      kdClusterMap[kdhit] = hit;
      kdClusters.push_back(kdhit);

      // Debug histogramming
      if (m_debugPlots && m_eventNumber == 0) {
        m_conformalEvents->Fill(kdhit->getU(), kdhit->getV());
        m_nonconformalEvents->Fill(hit->getPosition()[0], hit->getPosition()[1]);
        m_conformalEventsRTheta->Fill(kdhit->getR(), kdhit->getTheta());
      }
    }

    // Sort the KDClusters from larger to smaller radius
    streamlog_out(DEBUG4) << "Number of hits: " << kdClusters.size() << std::endl;
    std::sort(kdClusters.begin(), kdClusters.end(), sort_by_radiusKD);

    // Make the binary search tree. This tree class contains two binary trees - one sorted by u-v and the other by theta
    KDTree* nearestNeighbours = new KDTree(kdClusters);

    // Try to extend current tracks (if any) with unused hits
    int nCurrentTracks = finalCAtracks.size();
    streamlog_out(DEBUG4) << "Seeding with tracks" << std::endl;
    streamlog_out(DEBUG4) << "Attempting to extend current tracks: " << nCurrentTracks << std::endl;

    // Loop over all current tracks
    for (int currentTrack = 0; currentTrack < nCurrentTracks; currentTrack++) {
      // This step (although first) only runs when tracks have already been produced. An attempt
      // is made to extend them with hits from the new collection, using the final cell of the track
      // as the seed cell.

      // Containers to hold new cells made, and to check if a hit already has a cell connected to it
      std::vector<Cell*> cells;

      // Add the seed cell (the last cell added to the track) TODO: are we only creating mini-tracks and throwing away the first part of the track found? I think so!
      cells.push_back(finalCAtracks[currentTrack].front());

      // All seed cells have been created, now try create all "downstream" cells until no more can be added
      extendSeedCells(cells, used, nearestNeighbours, false);

      // We create all acceptable tracks by looping over all cells with high enough weight to create
      // a track and trace their route back to the seed hit. We then have to choose the best candidate
      // at the end (by minimum chi2 of a linear fit)
      std::map<Cell*, bool> usedCells;
      std::vector<cellularTrack> cellTracks;

      // Sort Cells from highest to lowest weight
      std::sort(cells.begin(), cells.end(), sort_by_cellWeight);

      // Create tracks by following a path along cells
      int nCells = cells.size();
      for (int itCell = 0; itCell < nCells; itCell++) {
        // Check if this cell has already been used
        if (usedCells.count(cells[itCell]))
          continue;

        // Produce all tracks leading back to the seed hit from this cell
        std::vector<cellularTrack> candidateTracks = createTracks(cells[itCell], usedCells);

        if (candidateTracks.size() == 0)
          continue;
        cellularTrack bestTrack;
        if (candidateTracks.size() == 1) {
          bestTrack = candidateTracks[0];
        } else {
          bestTrack = getLowestChi2(candidateTracks);
        }

        // Store track for later TODO: is this necessary? Twice performing best chi2... (first for all tracks from one seed, second for all "best" tracks)
        cellTracks.push_back(bestTrack);

        for (unsigned int trackCell = 0; trackCell < bestTrack.size(); trackCell++) {
          usedCells[bestTrack[trackCell]] = true;
        }
      }

      // CUT ON CHI2/NDOF of straight line fit here
      if (cellTracks.size() == 0)
        continue;
      cellularTrack bestTrack;
      if (cellTracks.size() == 1) {
        bestTrack = cellTracks[0];
      } else {
        bestTrack = getLowestChi2(cellTracks);
      }

      // Mark all hits as having been used
      KDCluster* kdStart = bestTrack[0]->getStart();
      used[kdStart]      = true;
      for (unsigned int trackCell = 0; trackCell < bestTrack.size(); trackCell++) {
        KDCluster* kdEnd = bestTrack[trackCell]->getEnd();
        used[kdEnd]      = true;
      }

      // Store the track
      finalCAtracks[currentTrack].clear();
      finalCAtracks[currentTrack].insert(finalCAtracks[currentTrack].end(), bestTrack.begin(), bestTrack.end());
    }

    // Try to create new tracks using all of the kdHits currently held
    unsigned int nKDHits = kdClusters.size();
    streamlog_out(DEBUG4) << "Seeding with hits" << std::endl;
    streamlog_out(DEBUG4) << "Attempting to seed with hits: " << nKDHits << std::endl;

    // Loop over all current hits
    for (unsigned int nKDHit = 0; nKDHit < nKDHits; nKDHit++) {
      // Get the kdHit and check if it has already been used (assigned to a track)
      KDCluster* kdhit = kdClusters[nKDHit];
      if (used.count(kdhit))
        continue;
      if (kdhit->getR() < 0.001)
        break;  // new cut - once we get to inner radius we will never make tracks. temp? TODO: make parameter?

      // The tracking differentiates between the first and all subsequent hits on a chain.
      // First, take the seed hit and look for sensible hits nearby to make an initial set
      // of cells. Once these are found, extrapolate the cells and look for additional hits
      // to produce a new cell along the chain. This is done mainly for speed: you ignore
      // all combinations which would be excluded when compared with the seed cell. In the
      // end we will try to produce a track starting from this seed, then begin again for
      // the next seed hit.

      // Get the initial seed cells for this hit, by looking at neighbouring hits in a given
      // theta slice around the hit. Check that they are at lower radius and further from the IP
      double     theta = kdhit->getTheta();
      VecCluster results;
      //      nearestNeighbours->allNeighboursInRadius(kdhit, m_maxDistance, results);
      nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);

      // Sort the neighbours from outer to inner radius
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
        if (used.count(nhit))
          continue;
        if (kdhit->sameLayer(nhit))
          continue;
        if (nhit->getR() >= kdhit->getR())
          continue;
        if (kdhit->getZ() > 0 && nhit->getZ() < kdhit->getZ())
          continue;
        if (kdhit->getZ() < 0 && nhit->getZ() > kdhit->getZ())
          continue;

        // Once we are far enough away we can never fulfil the distance requirements
        if ((kdhit->getR() - nhit->getR()) > m_maxDistance)
          break;

        // If nearest neighbours in theta, must put a distance cut
        double distance2 = (nhit->getU() - kdhit->getU()) * (nhit->getU() - kdhit->getU()) +
                           (nhit->getV() - kdhit->getV()) * (nhit->getV() - kdhit->getV());
        if (distance2 > (m_maxDistance * m_maxDistance))
          continue;

        // Create the new seed cell
        Cell* cell = new Cell(kdhit, nhit);
        cells.push_back(cell);
      }

      // All seed cells have been created, now try create all "downstream" cells until no more can be added
      extendSeedCells(cells, used, nearestNeighbours, false);

      // Now have all cells stemming from this seed hit. If it is possible to produce a track (ie. cells with depth X) then we will now...
      //      if(depth < (m_minClustersOnTrack-1)) continue; // TODO: check if this is correct

      // We create all acceptable tracks by looping over all cells with high enough weight to create
      // a track and trace their route back to the seed hit. We then have to choose the best candidate
      // at the end (by minimum chi2 of a linear fit)
      std::map<Cell*, bool> usedCells;
      std::vector<cellularTrack> cellTracks;

      // Sort Cells from highest to lowest weight
      std::sort(cells.begin(), cells.end(), sort_by_cellWeight);

      // Create tracks by following a path along cells
      int nCells = cells.size();
      for (int itCell = 0; itCell < nCells; itCell++) {
        // Check if this cell could produce a track (is on a long enough chain)
        if (cells[itCell]->getWeight() < (m_minClustersOnTrack - 2))
          break;

        // Check if this cell has already been used
        if (usedCells.count(cells[itCell]))
          continue;

        // Produce all tracks leading back to the seed hit from this cell
        std::vector<cellularTrack> candidateTracks = createTracks(cells[itCell], usedCells);

        if (candidateTracks.size() == 0)
          continue;
        cellularTrack bestTrack;
        if (candidateTracks.size() == 1) {
          bestTrack = candidateTracks[0];
        } else {
          bestTrack = getLowestChi2(candidateTracks);
        }

        // Store track for later TODO: is this necessary? Twice performing best chi2... (first for all tracks from one seed, second for all "best" tracks)
        cellTracks.push_back(bestTrack);

        for (unsigned int trackCell = 0; trackCell < bestTrack.size(); trackCell++) {
          usedCells[bestTrack[trackCell]] = true;
        }
      }

      // CUT ON CHI2/NDOF of straight line fit here
      if (cellTracks.size() == 0)
        continue;
      cellularTrack bestTrack;
      if (cellTracks.size() == 1) {
        bestTrack = cellTracks[0];
      } else {
        bestTrack = getLowestChi2(cellTracks);
      }

      // Mark all hits as having been used
      KDCluster* kdStart = bestTrack[0]->getStart();
      used[kdStart]      = true;
      for (unsigned int trackCell = 0; trackCell < bestTrack.size(); trackCell++) {
        KDCluster* kdEnd = bestTrack[trackCell]->getEnd();
        used[kdEnd]      = true;
      }

      // Store the track
      finalCAtracks.push_back(bestTrack);
    }

    // Now finished looking at this collection. All tracks which can have extra hits added
    // have had them added, and no new tracks were possible using the sum of all collections
    // till now. Add the next collection and try again...
    delete nearestNeighbours;
  }

  // Now make tracks from all of the candidates
  std::cout << "*** CA has made " << finalCAtracks.size() << " tracks ***" << std::endl;

  // Loop over all track candidates
  for (unsigned int caTrack = 0; caTrack < finalCAtracks.size(); caTrack++) {
    // Vector of all the hits on the track
    std::vector<TrackerHit*> trackHits;

    // Each CA track contains a list of cells. Get the kdHits, and the TrackerHits they point to
    KDCluster* kdEnd = finalCAtracks[caTrack][0]->getEnd();
    trackHits.push_back(kdClusterMap[kdEnd]);

    // Debug info on hit added
    streamlog_out(DEBUG5) << "  -- cluster 0: " << kdEnd->getU() << "," << kdEnd->getV() << " from detector "
                          << kdEnd->getSubdetector() << std::endl;

    // Loop over all cells and get the hit that they connect to
    for (unsigned int trackCell = 0; trackCell < finalCAtracks[caTrack].size(); trackCell++) {
      KDCluster* kdStart = finalCAtracks[caTrack][trackCell]->getStart();

      // Debug info on hit added
      streamlog_out(DEBUG5) << "  -- cluster " << trackCell + 1 << ": " << kdStart->getU() << "," << kdStart->getV()
                            << " from detector " << kdStart->getSubdetector() << std::endl;

      // Store the track hit that corresponds to this kdHit
      trackHits.push_back(kdClusterMap[kdStart]);

      // Debug plotting
      if (m_debugPlots) {
        // Draw the cells for event 0
        if (m_eventNumber == 0) {
          m_canvConformalEventDisplay->cd();
          drawline(finalCAtracks[caTrack][trackCell]->getStart(), finalCAtracks[caTrack][trackCell]->getEnd(),
                   finalCAtracks[caTrack].size() - trackCell);
        }

        // Fill the debug/tuning plots
        if (trackCell == finalCAtracks[caTrack].size() - 1)
          continue;
        double     angleBetweenCells = finalCAtracks[caTrack][trackCell]->getAngle(finalCAtracks[caTrack][trackCell + 1]);
        KDCluster* cluster0          = finalCAtracks[caTrack][trackCell]->getStart();
        KDCluster* cluster1          = finalCAtracks[caTrack][trackCell]->getEnd();
        KDCluster* cluster2          = finalCAtracks[caTrack][trackCell + 1]->getEnd();

        double cell0Length = sqrt(pow(cluster0->getU() - cluster1->getU(), 2) + pow(cluster0->getV() - cluster1->getV(), 2));
        double cell1Length = sqrt(pow(cluster1->getU() - cluster2->getU(), 2) + pow(cluster1->getV() - cluster2->getV(), 2));

        m_cellAngle->Fill(angleBetweenCells);
        m_cellAngleRadius->Fill(cluster2->getR(), angleBetweenCells);
        m_cellLengthRadius->Fill(cluster0->getR(), cell0Length);
        m_cellAngleLength->Fill(cell1Length, angleBetweenCells);
      }
    }

    // Only make tracks with n or more hits
    if (trackHits.size() < (unsigned int)m_minClustersOnTrack)
      continue;
    streamlog_out(DEBUG5) << "Made a track with " << trackHits.size() << " hits" << std::endl;

    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(), trackHits.end(), sort_by_radius);

    // Now we can make the track object and relations object, and fit the track
    TrackImpl* track = new TrackImpl;

    // First, for some reason there are 2 track objects, one which gets saved and one which is used for fitting. Don't ask...
    MarlinTrk::IMarlinTrack* marlinTrack = trackFactory->createTrack();

    // Save a vector of the hits to be used (why is this not attached to the track directly?? MarlinTrkUtils to be updated?)
    EVENT::TrackerHitVec trackfitHits;
    for (unsigned int itTrackHit = 0; itTrackHit < trackHits.size(); itTrackHit++) {
      trackfitHits.push_back(trackHits[itTrackHit]);
    }

    // Make an initial covariance matrix with very broad default values
    EVENT::FloatVec covMatrix(15, 0);             // Size 15, filled with 0s
    covMatrix[0]  = (m_initialTrackError_d0);     //sigma_d0^2
    covMatrix[2]  = (m_initialTrackError_phi0);   //sigma_phi0^2
    covMatrix[5]  = (m_initialTrackError_omega);  //sigma_omega^2
    covMatrix[9]  = (m_initialTrackError_z0);     //sigma_z0^2
    covMatrix[14] = (m_initialTrackError_tanL);   //sigma_tanl^2

    // Try to fit
    int fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackfitHits, track, MarlinTrk::IMarlinTrack::forward,
                                                       covMatrix, m_magneticField, m_maxChi2perHit);

    // Check track quality - if fit fails chi2 will be 0. For the moment add hits by hand to any track that fails the track fit, and store it as if it were ok...
    if (track->getChi2() == 0.) {
      std::cout << "Fit failed. Track has " << track->getTrackerHits().size() << " hits" << std::endl;
      std::cout << "Fit fail error " << fitError << std::endl;
      for (unsigned int p = 0; p < trackfitHits.size(); p++) {
        track->addHit(trackfitHits[p]);
      }
    }  // delete track; delete marlinTrack; continue;}

    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2] = trackHits.size();

    // Push back to the output container
    trackCollection->addElement(track);
  }

  // Draw the conformal event display hits for debugging
  if (m_debugPlots && m_eventNumber == 0) {
    m_canvConformalEventDisplay->cd();
    m_conformalEvents->DrawCopy("same");
  }

  // Save the output track collection
  evt->addCollection(trackCollection, m_outputTrackCollection);

  // Increment the event number
  m_eventNumber++;
}

void ConformalTracking::check(LCEvent* evt) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void ConformalTracking::end() {
  streamlog_out(MESSAGE) << " end()  " << name() << " processed " << m_eventNumber << " events in " << m_runNumber
                         << " runs " << std::endl;

  m_canvConformalEventDisplay->Write();
  m_canvConformalEventDisplayMC->Write();
}

void ConformalTracking::getCollection(LCCollection*& collection, std::string collectionName, LCEvent* evt) {
  try {
    collection = evt->getCollection(collectionName);
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void ConformalTracking::extendSeedCells(std::vector<Cell*>& cells, std::map<KDCluster*, bool> used,
                                        KDTree* nearestNeighbours, bool extendingTrack) {
  unsigned int nCells   = 0;
  int          depth    = 0;
  int          startPos = 0;

  // If extending an existing track, don't redo the search from the beginning of the track. Just extend from
  // the last cell
  if (extendingTrack) {
    depth    = cells.size() - 1;
    startPos = cells.size() - 1;
  }

  // Keep track of existing cells in case there are branches in the track
  std::map<KDCluster*, std::vector<Cell*>> existingCells;

  // Try to create all "downstream" cells until no more can be added
  while (cells.size() != nCells) {
    // Extend all cells with depth N. In the next iteration, look at cells with depth N+1
    nCells = cells.size();
    for (unsigned int itCell = startPos; itCell < nCells; itCell++) {
      // Get the end point of the cell (to search for neighbouring hits to form new cells connected to this one)
      KDCluster* hit = cells[itCell]->getEnd();

      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      KDCluster* fakeHit =
          extrapolateCell(cells[itCell], m_maxDistance / 2.);  // TODO: make this search a function of radius
      VecCluster results;
      nearestNeighbours->allNeighboursInRadius(fakeHit, m_maxDistance / 2., results);
      delete fakeHit;

      // Make new cells pointing inwards
      for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
        // Get the neighbouring hit
        KDCluster* nhit = results[neighbour];

        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
        if (used.count(nhit))
          continue;
        if (hit->sameLayer(nhit))
          continue;
        if (nhit->getR() >= hit->getR())
          continue;
        if (hit->getZ() > 0 && nhit->getZ() < hit->getZ())
          continue;
        if (hit->getZ() < 0 && nhit->getZ() > hit->getZ())
          continue;

        // Check if this cell already exists (rejoining branch)
        if (existingCells.count(hit) != 0) {
          bool alreadyExists  = false;
          int  nExistingCells = existingCells[hit].size();
          for (int iterExisting = 0; iterExisting < nExistingCells; iterExisting++) {
            if (existingCells[hit][iterExisting]->getEnd() == nhit) {
              alreadyExists = true;
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
        if (cells[itCell]->getAngle(cell) > (m_maxCellAngle * exp(-0.0015 / nhit->getR()))) {
          delete cell;
          continue;
        }

        // Set the information about which cell this new cell is attached to and store it
        cell->setFrom(cells[itCell]);
        cells[itCell]->setTo(cell);
        cells.push_back(cell);
        existingCells[hit].push_back(cell);
      }

      // Finished adding new cells to this cell
    }

    // All new cells added at this depth
    startPos = nCells;
    depth++;
  }

  // No more downstream cells can be added
}

void ConformalTracking::drawline(KDCluster* hitStart, KDCluster* hitEnd, int colour) {
  //  std::cout<<"Plotting initial cell"<<std::endl;

  TLine* line = new TLine(hitStart->getU(), hitStart->getV(), hitEnd->getU(), hitEnd->getV());
  line->SetLineColor(colour);
  line->Draw();
}

// TODO: check the logic in used cell notation, branching, weights of cells being followed. In particular if we increase weight and make continuously connected
// cells with non-continuous weight...
std::vector<cellularTrack> ConformalTracking::createTracks(Cell* seedCell, std::map<Cell*, bool>& usedCells) {
  //  std::cout<<" - Creating tracks"<<std::endl;
  cellularTrack seedTrack;
  seedTrack.push_back(seedCell);

  std::vector<cellularTrack> cellularTracks;
  cellularTracks.push_back(seedTrack);

  while (toBeUpdated(cellularTracks)) {
    int removeTrack = -1;
    int nTracks     = cellularTracks.size();
    //    std::cout<<" -- check for update on "<<nTracks<<" tracks"<<std::endl;
    for (int iTrack = 0; iTrack < nTracks; iTrack++) {
      if (removeTrack > (-1))
        continue;
      if (cellularTracks[iTrack].back()->getWeight() > 0)
        followPath(cellularTracks, iTrack, usedCells, removeTrack);
    }
    if (removeTrack > (-1)) {
      //      std::cout<<"Erasing track "<<removeTrack<<std::endl;
      cellularTracks.erase(cellularTracks.begin() + removeTrack);
    }
  }

  return cellularTracks;
}

bool ConformalTracking::toBeUpdated(std::vector<cellularTrack> cellularTracks) {
  bool update = false;
  for (unsigned int iTrack = 0; iTrack < cellularTracks.size(); iTrack++)
    if (cellularTracks[iTrack].back()->getWeight() > 0) {
      update = true;
      break;
    }
  return update;
}

void ConformalTracking::followPath(std::vector<cellularTrack>& cellularTracks, int trackNumber,
                                   std::map<Cell*, bool>& usedCells, int& removeTrack) {
  Cell* cell = cellularTracks[trackNumber].back();

  if (cell->getWeight() == 0)
    return;

  //  std::cout<<"  -- following path of cell with weight "<<cell->getWeight()<<std::endl;
  //  std::cout<<"  -- this cell is connected to "<<cell->getFrom().size()<<std::endl;
  if (cell->getFrom().size() > 1) {
    //    std::cout<<"  -- branching!"<<std::endl;
    if (usedCells.count(cell->getFrom()[0])) {
      removeTrack = trackNumber;
    } else {
      if ((cell->getWeight() - cell->getFrom()[0]->getWeight()) != 1) {
        removeTrack = trackNumber;
      } else {
        cellularTracks[trackNumber].push_back(cell->getFrom()[0]);
        //        std::cout<<"  -- created branch 0"<<std::endl;
      }
    }
    for (unsigned int itCell = 1; itCell < cell->getFrom().size(); itCell++) {
      if (usedCells.count(cell->getFrom()[itCell])) {
        continue;
      } else {
        if ((cell->getWeight() - cell->getFrom()[itCell]->getWeight()) != 1)
          continue;
        //        std::cout<<"  -- created branch "<<itCell<<std::endl;
        //        std::cout<<"  -- cell weight is "<<cell->getWeight()<<" and branch has weight "<<cell->getFrom()[itCell]->getWeight()<<std::endl;
        cellularTrack branchedTrack = cellularTracks[trackNumber];
        branchedTrack.push_back(cell->getFrom()[itCell]);
        cellularTracks.push_back(branchedTrack);
      }
    }
  }

  while (cell->getWeight() > 0 && cell->getFrom().size() == 1) {
    //    std::cout<<"  -- following linear path. Current weight "<<cell->getWeight()<<std::endl;
    Cell* parentCell = cell->getFrom()[0];
    if (usedCells.count(cell)) {
      removeTrack = trackNumber;
      break;
    }
    cellularTracks[trackNumber].push_back(parentCell);
    cell = parentCell;
  }

  //	if(cell->getFrom().size() == 0) return;

  return;
}

cellularTrack ConformalTracking::getLowestChi2(std::vector<cellularTrack> candidateTracks) {
  //  std::cout<<"GETTING lowest chi2"<<std::endl;
  double bestChi2ndof = 1000000.;
  int    bestTrack    = -1;
  for (unsigned int itTrack = 0; itTrack < candidateTracks.size(); itTrack++) {
    //    std::cout<<" - track "<<itTrack<<std::endl;
    TLinearFitter* fitter = new TLinearFitter(1);

    fitter->SetFormula("pol1");
    double     npoints = 0.;
    KDCluster* kdStart = candidateTracks[itTrack][0]->getStart();
    double*    u       = new double[0];
    u[0]               = kdStart->getU();
    fitter->AddPoint(u, kdStart->getV(), 1. / kdStart->getR());
    npoints++;

    for (unsigned int trackCell = 0; trackCell < candidateTracks[itTrack].size(); trackCell++) {
      KDCluster* kdEnd = candidateTracks[itTrack][trackCell]->getEnd();
      u[0]             = kdEnd->getU();
      fitter->AddPoint(u, kdEnd->getV(), 1. / kdEnd->getR());
      npoints++;
    }

    fitter->Eval();
    double chi2ndof = fitter->GetChisquare() / (npoints - 2.);
    //    std::cout<<" - chi2/ndof = "<<chi2ndof<<std::endl;
    if (chi2ndof < bestChi2ndof) {
      bestChi2ndof = chi2ndof;
      bestTrack    = itTrack;
    }
    delete fitter;
  }
  //  std::cout<<" - returning track "<<bestTrack<<std::endl;
  return candidateTracks[bestTrack];
}

void ConformalTracking::updateCell(Cell* cell) {
  if (cell->getTo().size() != 0) {
    for (unsigned int i = 0; i < cell->getTo().size(); i++) {
      cell->getTo()[i]->update(cell);
      updateCell(cell->getTo()[i]);
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
