#include "ConformalTracking.h"
#include "DDRec/API/IDDecoder.h"

#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>

#include <UTIL/CellIDEncoder.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"

#include "CLHEP/Vector/TwoVector.h"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include "TLinearFitter.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"
#include "TLinearFitter.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"

#include <cmath>
#include <map>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>

#include "KDTree.h"
#include "Cell.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace DD4hep ;
using namespace AIDA ;

ConformalTracking aConformalTracking;

/*
 
 Pattern recognition code for the CLIC detector, using conformal mapping and cellular automaton
 
 */

ConformalTracking::ConformalTracking() : Processor("ConformalTracking") {
  
  // Processor description
  _description = "ConformalTracking constructs tracks using a combined conformal mapping and cellular automaton approach." ;
  
  // Input collections - tracker hits
  std::vector<std::string> inputTrackerHitCollections;
  inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("VXDEndcapTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerEndcapHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerEndcapHits"));
  registerInputCollections( LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames" ,"Name of the TrackerHit input collections", m_inputTrackerHitCollections, inputTrackerHitCollections ) ;
  
  // Debugging collections - MC particles and relation collections
  registerInputCollection( LCIO::MCPARTICLE, "MCParticleCollectionName", "Name of the MCParticle input collection", m_inputParticleCollection, std::string("MCParticle"));
  std::vector<std::string> inputRelationCollections;
  inputRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
  inputRelationCollections.push_back(std::string("VXDEndcapTrackerHitRelations"));
  inputRelationCollections.push_back(std::string("InnerTrackerBarrelHitsRelations"));
  inputRelationCollections.push_back(std::string("OuterTrackerBarrelHitsRelations"));
  inputRelationCollections.push_back(std::string("InnerTrackerEndcapHitsRelations"));
  inputRelationCollections.push_back(std::string("OuterTrackerEndcapHitsRelations"));
  registerInputCollections( LCIO::LCRELATION, "RelationsNames" ,"Name of the TrackerHit relation collections", m_inputRelationCollections, inputRelationCollections ) ;
  
  // Output collections - tracks
  registerOutputCollection( LCIO::TRACK, "SiTrackCollectionName", "Silicon track Collection Name", m_outputTrackCollection, std::string("CATracks"));
  registerOutputCollection( LCIO::TRACKERHITPLANE, "DebugHits", "DebugHits", m_outputDebugHits, std::string("DebugHits"));
  
  // Parameters for tracking
  registerProcessorParameter( "DebugPlots",         "Plots for debugging the tracking",				  m_debugPlots,		  bool(false) 	);
  registerProcessorParameter( "ThetaRange", 	    "Angular range for initial cell seeding",			  m_thetaRange,		  double(0.1)	);
  registerProcessorParameter( "MaxCellAngle", 	    "Cut on angle between two cells for cell to be valid",	  m_maxCellAngle,	  double(0.035)	);
  registerProcessorParameter( "MaxCellAngleRZ",     "Cut on angle between two cells in RZ for cell to be valid",  m_maxCellAngleRZ, 	  double(0.035)	);
  registerProcessorParameter( "MaxDistance",        "Maximum length of a cell (max. distance between two hits)",  m_maxDistance,	  double(0.015) );
  registerProcessorParameter( "MaxChi2",            "Maximum chi2/ndof for linear conformal tracks",		  m_chi2cut,		  double(300.) 	);
  registerProcessorParameter( "MinClustersOnTrack", "Minimum number of clusters to create a track",		  m_minClustersOnTrack,	  int(6)	);
  registerProcessorParameter( "trackPurity",        "Purity value used for checking if tracks are real or not",	  m_purity,               double(0.75)	);

}

// Sort tracker hits from smaller to larger radius
bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2){
  double radius1 = sqrt((hit1->getPosition()[0])*(hit1->getPosition()[0]) + (hit1->getPosition()[1])*(hit1->getPosition()[1]));
  double radius2 = sqrt((hit2->getPosition()[0])*(hit2->getPosition()[0]) + (hit2->getPosition()[1])*(hit2->getPosition()[1]));
  return (radius1 < radius2);
}

// Sort kd hits from larger to smaller radius
bool sort_by_radiusKD(KDCluster* hit1, KDCluster* hit2){
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 > radius2);
}

// Sort kd hits from smaller to larger radius
bool sort_by_lower_radiusKD(KDCluster* hit1, KDCluster* hit2){
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 < radius2);
}

// Sort kdhits by lower to higher layer number
bool sort_by_layer(KDCluster* hit1, KDCluster* hit2){
  
  if( hit1->getSubdetector() != hit2->getSubdetector() ) return (hit1->getSubdetector() < hit2->getSubdetector());
  else if( hit1->getSide() != hit2->getSide() ) return (hit1->getSide() < hit2->getSide());
  else if( hit1->getLayer() != hit2->getLayer() ) return (hit1->getLayer() < hit2->getLayer());
  else return false;
}

// Sort cells from higher to lower weight
bool sort_by_cellWeight(Cell* cell1, Cell* cell2){
  int weight1 = cell1->getWeight();
  int weight2 = cell2->getWeight();
  return (weight1 > weight2);
}

// Sort kdtracks from longest to shortest
bool sort_by_length(KDTrack* track1, KDTrack* track2){
  return (track1->m_clusters.size() > track2->m_clusters.size());
}


void ConformalTracking::init() {
  
  // Print the initial parameters
  printParameters() ;
  
  // Reset counters
  m_runNumber = 0 ;
  m_eventNumber = 0 ;
  
  // Set up the track fit factory
  const gear::GearMgr* fakeGear = 0;
  trackFactory =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , fakeGear , "" ) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        true) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       true) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  false) ;
  trackFactory->init() ;
  
  // Put default values for track fitting
  m_initialTrackError_d0 = 1.e6;
  m_initialTrackError_phi0 = 1.e2;
  m_initialTrackError_omega = 1.e-4;
  m_initialTrackError_z0 = 1.e6;
  m_initialTrackError_tanL = 1.e2;
  m_maxChi2perHit = 1.e2;
  
  // Get the magnetic field
  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
  m_magneticField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  
  // Seed hit for debug printouts. If not set later, isn't used
  debugSeed = NULL;

  // Initialise histograms (if debug plotting on)
  if(m_debugPlots){
    
    // Automatically add histograms to output root file
    AIDAProcessor::histogramFactory(this);
    
    // Histogram initailisation
    m_szDistribution = new TH2F("m_szDistribution","m_szDistribution",20000,-100,100,200,-10,10);
    m_uvDistribution = new TH2F("m_uvDistribution","m_uvDistribution",1000,-0.05,0.05,1000,-0.05,0.05);
    m_xyDistribution = new TH2F("m_xyDistribution","m_xyDistribution",500,-1500,1500,500,-1500,1500);
    m_xyzDistribution = new TH3F("m_xyzDistribution","m_xyzDistribution",50,0,100,50,0,100,100,0,25);

    // Histograms for tuning parameters (cell angle cut, cell length cut)
    m_cellAngle = new TH1F("cellAngle","cellAngle",1250,0,0.05);
    m_cellAngleRadius = new TH2F("cellAngleRadius","cellAngleRadius",400,0,0.04,1000,0,0.04);
    m_cellLengthRadius = new TH2F("cellLengthRadius","cellLengthRadius",300,0,0.03,1000,0,0.04);
    m_cellAngleLength = new TH2F("cellAngleLength","cellAngleLength",400,0,0.04,300,0,0.03);
    m_conformalChi2 = new TH1F("conformalChi2","conformalChi2",100,0,100);
    m_conformalChi2real = new TH1F("conformalChi2real","conformalChi2real",1000,0,1000);
    m_conformalChi2fake = new TH1F("conformalChi2fake","conformalChi2fake",1000,0,1000);
    m_conformalChi2Purity = new TH2F("conformalChi2Purity","conformalChi2Purity",150,0,1.5,1000,0,1000);
    m_conformalChi2MC = new TH1F("conformalChi2MC","conformalChi2MC",1000,0,1000);
    
    m_cellAngleMC = new TH1F("cellAngleMC","cellAngleMC",1250,0,0.05);
    m_cellAngleRadiusMC = new TH2F("cellAngleRadiusMC","cellAngleRadiusMC",400,0,0.04,1000,0,0.04);
    m_cellLengthRadiusMC = new TH2F("cellLengthRadiusMC","cellLengthRadiusMC",300,0,0.03,1000,0,0.04);
    m_cellAngleLengthMC = new TH2F("cellAngleLengthMC","cellAngleLengthMC",400,0,0.04,300,0,0.03);

    m_cellAngleRZMC = new TH1F("cellAngleRZMC","cellAngleRZMC",1250,0,0.05);

    // Histograms for "event display"
    m_conformalEvents = new TH2F("conformalEvents","conformalEvents",1000,-0.05,0.05,1000,-0.05,0.05);
    m_nonconformalEvents = new TH2F("nonconformalEvents","nonconformalEvents",500,-1500,1500,500,-1500,1500);
    m_conformalEventsRTheta = new TH2F("conformalEventsRTheta","conformalEventsRTheta",200,0,0.05,632,-0.02,6.30);
    m_conformalEventsMC = new TH2F("conformalEventsMC","conformalEventsMC",1000,-0.05,0.05,1000,-0.05,0.05);
    
    m_canvConformalEventDisplay = new TCanvas("canvConformalEventDisplay","canvConformalEventDisplay");
    m_canvConformalEventDisplayAllCells = new TCanvas("canvConformalEventDisplayAllCells","canvConformalEventDisplayAllCells");
    m_canvConformalEventDisplayAcceptedCells = new TCanvas("canvConformalEventDisplayAcceptedCells","canvConformalEventDisplayAcceptedCells");
    m_canvConformalEventDisplayMC = new TCanvas("canvConformalEventDisplayMC","canvConformalEventDisplayMC");
    m_canvConformalEventDisplayMCunreconstructed = new TCanvas("canvConformalEventDisplayMCunreconstructed","canvConformalEventDisplayMCunreconstructed");
    
  }
  
  // Register this process
  Global::EVENTSEEDER->registerProcessor(this);
		
}


void ConformalTracking::processRunHeader( LCRunHeader* run) {
  ++m_runNumber ;
}

// The main code, run over each event
void ConformalTracking::processEvent( LCEvent* evt ) {
  
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
  
  streamlog_out( DEBUG4 )<<"Event number: "<<m_eventNumber<<std::endl;
  
  // Object to store all of the track hit collections passed to the pattern recognition
  std::vector<LCCollection*> trackerHitCollections;
  std::vector<LCRelationNavigator*> relations;

  // Loop over each input collection and get the hits
  for(unsigned int collection=0; collection<m_inputTrackerHitCollections.size();collection++){
    
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = 0 ;
    getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt); if(trackerHitCollection == 0) continue;
    streamlog_out( DEBUG4 )<<"Collection "<<m_inputTrackerHitCollections[collection]<<" contains "<<trackerHitCollection->getNumberOfElements()<<" hits"<<std::endl;
    trackerHitCollections.push_back(trackerHitCollection);
    
    // If debugging, get the relations between tracker hits and MC particle
    if(m_debugPlots){
      // Get the collection of tracker hit relations
      LCCollection* trackerHitRelationCollection = 0 ;
      getCollection(trackerHitRelationCollection, m_inputRelationCollections[collection], evt); if(trackerHitRelationCollection == 0) continue;
      // Create the relations navigator
      LCRelationNavigator* relation = new LCRelationNavigator( trackerHitRelationCollection );
      relations.push_back(relation);
    }

  }
  
  // Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec( LCIO::TRACK )  ;
  LCCollectionVec* debugHitCollection = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
  debugHitCollection->setSubset(true);
  
  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackCollection->setFlag( trkFlag.getFlag()  ) ;
  
  // Set up ID decoder
  UTIL::BitField64 m_encoder( lcio::LCTrackerCellID::encoding_string() ) ;

  /*
   Debug plotting. This section picks up tracks reconstructed using the cheated pattern recognition (TruthTrackFinder) and uses it to show
   values which are cut on during the tracking. This can be used to tune the cut ranges.
   */
  
  std::vector<KDCluster*> debugHits;
  map<KDCluster*, TrackerHitPlane*> tempHolder;
  map<KDCluster*,MCParticle*> kdParticles;
  map<TrackerHitPlane*, KDCluster*> conformalHits;
  map<MCParticle*,bool> reconstructed;
  
    // Container to store the hits
  std::map<MCParticle*, std::vector<KDCluster*> > particleHits;
  LCCollection* particleCollection = 0 ;
  
  if(m_debugPlots){
    
    // Get the MC particle collection
    getCollection(particleCollection, m_inputParticleCollection, evt);
    if(particleCollection == 0){
      delete trackCollection;
      delete debugHitCollection;
      return;
    }
    
    // Draw the empty event display onto the canvas, so that cells can be added sequentially
    if(m_eventNumber == 0){
      m_canvConformalEventDisplayMC->cd();
      m_conformalEventsMC->DrawCopy("");
      m_canvConformalEventDisplayMCunreconstructed->cd();
      m_conformalEventsMC->DrawCopy("");    }
    
    // Loop over all hits and assign them to a MC particle. Then loop over all particles and make diagnostic plots as if
    // it was the real tracking, highlighting any cells which would fail to be produced
    
    // Loop over all input collections
    for(unsigned int collection=0; collection<trackerHitCollections.size();collection++){
      // Loop over tracker hits
      int nHits = trackerHitCollections[collection]->getNumberOfElements();
      for(int itHit=0;itHit<nHits;itHit++){
        // Get the hit
        TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
        // Get the related simulated hit(s)
        const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects( hit );
        // Take the first hit only (TODO: this should be changed? Loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit?)
        SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
        // Get the particle belonging to that hit
        MCParticle* particle = simHit->getMCParticle();
        // Make the conformal hit, first get subdetector information to check if it is a barrel or endcap
        const int celId = hit->getCellID0() ;
        m_encoder.setValue(celId) ;
        int side = m_encoder[lcio::LCTrackerCellID::side()];
        bool isEndcap = false;
        if(side != ILDDetID::barrel) isEndcap = true;
        KDCluster* kdhit = new KDCluster(hit,isEndcap);
        // Push back the element into the container
        particleHits[particle].push_back(kdhit);
        tempHolder[kdhit] = hit;
        kdParticles[kdhit] = particle;
        conformalHits[hit] = kdhit;
      }
    }
    
    // Now loop over all MC particles and make the cells connecting hits
    int nParticles = particleCollection->getNumberOfElements();
    for(int itP=0;itP<nParticles;itP++){
      // Get the particle
      MCParticle* mcParticle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itP) ) ;
      // Get the vector of hits from the container
      if(particleHits.count(mcParticle) == 0) continue;
      std::vector<KDCluster*> trackHits = particleHits[mcParticle];
      // Only make tracks with n or more hits
      if(trackHits.size() < (unsigned int)m_minClustersOnTrack) continue;
      // Discard low momentum particles
      double particlePt = sqrt( mcParticle->getMomentum()[0]*mcParticle->getMomentum()[0] + mcParticle->getMomentum()[1]*mcParticle->getMomentum()[1] );
      // Cut on stable particles
      if(mcParticle->getGeneratorStatus() != 1) continue;
      // Sort the hits from larger to smaller radius
      std::sort(trackHits.begin(),trackHits.end(),sort_by_radiusKD);
      
      // Make a track
      KDTrack* mcTrack = new KDTrack();
      // Loop over all hits for debugging
      for(int itHit=0;itHit<trackHits.size();itHit++){
        // Get the conformal clusters
        KDCluster* cluster = trackHits[itHit];
        mcTrack->add(cluster);
      }
      
      // Fit the track and plot the chi2
      mcTrack->linearRegression();
      m_conformalChi2MC->Fill(mcTrack->chi2ndof());
      delete mcTrack;
      
      // Now loop over the hits and make cells - filling histograms along the way
      int nHits = trackHits.size();
      for(int itHit=0;itHit<(nHits-2);itHit++){
        
        // Get the conformal clusters
        KDCluster* cluster0 = trackHits[itHit];
        KDCluster* cluster1 = trackHits[itHit+1];
        KDCluster* cluster2 = trackHits[itHit+2];
        
        // Make the two cells connecting these three hits
        Cell* cell = new Cell(cluster0,cluster1);
        cell->setWeight(itHit);
        Cell* cell1 = new Cell(cluster1,cluster2);
        cell1->setWeight(itHit+1);
        
        // Fill the debug/tuning plots
        double angleBetweenCells = cell->getAngle(cell1);
        double angleRZBetweenCells = cell->getAngleRZ(cell1);
        double cell0Length = sqrt( pow(cluster0->getU() - cluster1->getU(),2) + pow(cluster0->getV() - cluster1->getV(),2) );
        double cell1Length = sqrt( pow(cluster1->getU() - cluster2->getU(),2) + pow(cluster1->getV() - cluster2->getV(),2) );
        
        m_cellAngleMC->Fill(angleBetweenCells);
        m_cellAngleRadiusMC->Fill(cluster2->getR(),angleBetweenCells);
        m_cellLengthRadiusMC->Fill(cluster0->getR(),cell0Length);
        m_cellAngleLengthMC->Fill(cell1Length,angleBetweenCells);
        m_cellAngleRZMC->Fill(angleRZBetweenCells);
        
        // Draw cells on the first event
        if(m_eventNumber == 0){
          // Fill the event display (hit positions)
          m_conformalEventsMC->Fill( cluster0->getU(),cluster0->getV() );
          m_conformalEventsMC->Fill( cluster1->getU(),cluster1->getV() );
          m_conformalEventsMC->Fill( cluster2->getU(),cluster2->getV() );
          // Draw the cell lines on the event display. Use the line style to show
          // if the cells would have been cut by some of the search criteria
          m_canvConformalEventDisplayMC->cd();
          if(itHit == 0){
            drawline(cluster0,cluster1,itHit+1);
          }
          // Draw line style differently if the cell angle was too large
          if(angleBetweenCells > (m_maxCellAngle)){
            drawline(cluster1,cluster2,itHit+2,3);
          }else{
            drawline(cluster1,cluster2,itHit+2);
          }
        }
        delete cell; delete cell1;
      }
    }
    
    // Draw the final set of conformal hits (on top of the cell lines)
    if(m_eventNumber == 0){
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
  }
  
  /*
   Now start doing things!
   */
  
  // Some global containers to be used throughout the tracking. A collection of conformal hits will be made, with a link
  // pointing back to the corresponding cluster. A record of all used hits will be kept.
  
  std::map< int,std::vector<KDCluster*> > collectionClusters;	// Conformal hits
  std::map<KDCluster*,TrackerHitPlane*> kdClusterMap;			// Their link to "real" hits
  std::map<KDCluster*,bool> used;								// Map of whether a hit has been included in a track or not
  std::map<KDCluster*,bool> used2;
  std::vector<KDTrack*> conformalTracks;						// KD tracks - each is a list of kd hits in the found tracks

  // Create the conformal hit collections for each tracker hit collection (and save the link)
  for(unsigned int collection=0; collection<trackerHitCollections.size();collection++){
    
    // Loop over tracker hits and make conformal hit collection
    std::vector<KDCluster*> tempClusters;
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    for(int itHit=0;itHit<nHits;itHit++){
      
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
      
      // Get subdetector information and check if the hit is in the barrel or endcaps
      const int celId = hit->getCellID0() ;
      m_encoder.setValue(celId) ;
      int subdet = m_encoder[lcio::LCTrackerCellID::subdet()];
      int side = m_encoder[lcio::LCTrackerCellID::side()];
      int layer = m_encoder[lcio::LCTrackerCellID::layer()];
      bool isEndcap = false;
      if(side != ILDDetID::barrel) isEndcap = true;
      
      // Make a new kd cluster (if debugging then these already exist - don't remake them)
      KDCluster* kdhit;
      if(m_debugPlots) kdhit = conformalHits[hit];
      if(!m_debugPlots) kdhit = new KDCluster(hit,isEndcap);

      // Set the subdetector information
      kdhit->setDetectorInfo(subdet,side,layer);
      
      // Store the link between the two
      kdClusterMap[kdhit] = hit;
      tempClusters.push_back(kdhit);
      
      // Debug histogramming
      if(m_debugPlots && m_eventNumber == 0){
        m_conformalEvents->Fill( kdhit->getU(),kdhit->getV() );
        m_nonconformalEvents->Fill( hit->getPosition()[0], hit->getPosition()[1] );
        m_conformalEventsRTheta->Fill( kdhit->getR(), kdhit->getTheta() );
      }
      
    }
    
    collectionClusters[collection] = tempClusters;
  }
  
  // Loop over all input collections. Tracking will be attempted on the collection, then hits from the next collection
  // will be added to the unused hits already there.
  int seedCollections = 0;
  for(unsigned int collection=0; collection<trackerHitCollections.size();collection++){
  
    // The set of conformal hits which will be considered in this iteration
    std::vector<KDCluster*> kdClusters;

    if(collection < seedCollections){
      std::vector<KDCluster*> clusters = collectionClusters[collection]; //this makes a copy FIX ME
      int nhits = clusters.size();
      for(int hit=0;hit<nhits;hit++){
        kdClusters.push_back(clusters[hit]);
      }
    }else{
      // Add hits from this and previous collections to the list
    	for(unsigned int col=0;col<=collection;col++){
        std::vector<KDCluster*> clusters = collectionClusters[col]; //this makes a copy FIX ME
        int nhits = clusters.size();
        for(int hit=0;hit<nhits;hit++){
//          if( clusters[hit]->used() ) continue;
          kdClusters.push_back(clusters[hit]);
        }
      }
    }
    
    // Sort the KDClusters from larger to smaller radius
    if(kdClusters.size() == 0) continue;
    std::cout<<"Number of hits: "<<kdClusters.size()<<std::endl;
    streamlog_out( DEBUG4 )<<"Number of hits: "<<kdClusters.size()<<std::endl;
    std::sort(kdClusters.begin(),kdClusters.end(),sort_by_radiusKD);
    
    // Make the binary search tree. This tree class contains two binary trees - one sorted by u-v and the other by theta
    KDTree* nearestNeighbours = new KDTree(kdClusters,m_thetaRange);
  
  	// ---------------------------------------------------------------------
    // Try to extend current tracks (if any) with unused hits
  	// ---------------------------------------------------------------------
    int nCurrentTracks = conformalTracks.size();
    streamlog_out( DEBUG4 )<<"Seeding with tracks"<<std::endl;
    streamlog_out( DEBUG4 )<<"Attempting to extend current tracks: "<<nCurrentTracks<<std::endl;
//    std::cout<<"Seeding with tracks"<<std::endl;
    
    // Loop over all current tracks
/*    for(int currentTrack=0;currentTrack<nCurrentTracks;currentTrack++){

//      continue;
//      cout<<"== Trying to extend track "<<currentTrack<<endl;
      // This step (although first) only runs when tracks have already been produced. An attempt
      // is made to extend them with hits from the new collection, using the final cell of the track
      // as the seed cell.
      
      // Containers to hold new cells made, and to check if a hit already has a cell connected to it
      std::vector<Cell*> cells;
      
      // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
      Cell* seedCell = new Cell(conformalTracks[currentTrack]->m_clusters[1],conformalTracks[currentTrack]->m_clusters[0]);
      cells.push_back(seedCell);
            
      // All seed cells have been created, now try create all "downstream" cells until no more can be added
      extendSeedCells(cells, used, nearestNeighbours, false, debugHits);
      
      // We create all acceptable tracks by looping over all cells with high enough weight to create
      // a track and trace their route back to the seed hit. We then have to choose the best candidate
      // at the end (by minimum chi2 of a linear fit)
      std::map<Cell*,bool> usedCells;
      std::map<Cell*,bool> usedCells2;
      std::vector<cellularTrack*> trackSegments;
      
      // Sort Cells from highest to lowest weight
      std::sort(cells.begin(),cells.end(),sort_by_cellWeight);
      
      // Create track "segments" leading back to the current track being considered
      int nCells = cells.size();
      for(int itCell=0;itCell<nCells;itCell++){

        // Check if this cell has already been used
        if(usedCells.count(cells[itCell])) continue;
        
        // Decide if we this cell has enough hits on it. In order to avoid spurious hit addition
        // we would like at least 2 new hits to be added to the track. This means a cell weight of
        // at least 2
        if(cells[itCell]->getWeight() < 1) break;
        
        // Produce all segments leading back to the track from this cell
        std::vector<cellularTrack*> candidateSegments;
        createTracksNew(candidateSegments,cells[itCell],usedCells2);
        
        // Store all of these segments for later
        if(candidateSegments.size() == 0) continue;
        trackSegments.insert(trackSegments.end(),candidateSegments.begin(), candidateSegments.end());

        // Mark the cells from these segments as having been used
//        for(unsigned int itSegment=0;itSegment<candidateSegments.size();itSegment++){
//          for(unsigned int itCell=0;itCell<candidateSegments[itSegment].size();itCell++){
//            usedCells[candidateSegments[itSegment][itCell]]=true;
//          }
//        }
      }

      // Decide which segment to add on to the track, and mark the added hits as used
      if(trackSegments.size() == 0) continue;
      extendTrack(conformalTracks[currentTrack],trackSegments,used,usedCells);
      
      // Clean up
      for(unsigned int itCell=0;itCell<cells.size();itCell++) delete cells[itCell];

    }
*/
    // Loop over all current tracks
    std::cout<<"EXTENDING tracks"<<std::endl;
    
    std::vector<KDCluster*> tempClusterContainer = kdClusters;
    std::sort(tempClusterContainer.begin(), tempClusterContainer.end(), sort_by_layer);

    for(int currentTrack=0;currentTrack<nCurrentTracks;currentTrack++){
      
      if(collection < 1) continue;
      
      // Make sure that all tracks have a kalman track attached to them
      if(conformalTracks[currentTrack]->kalmanTrack() == NULL){
        KalmanTrack* kalmanTrack = new KalmanTrack(conformalTracks[currentTrack]);
        conformalTracks[currentTrack]->setKalmanTrack(kalmanTrack);
      }
      
      // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
      Cell* seedCell = new Cell(conformalTracks[currentTrack]->m_clusters[1],conformalTracks[currentTrack]->m_clusters[0]);

      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      double searchDistance = m_maxDistance;
      KDCluster* fakeHit = extrapolateCell(seedCell,searchDistance/2.); // TODO: make this search a function of radius
      VecCluster results2;
      nearestNeighbours->allNeighboursInRadius(fakeHit, 10*searchDistance/2., results2);
      std::sort(results2.begin(),results2.end(),sort_by_radiusKD);
//      std::cout<<"- fake hit has position ("<<fakeHit->getU()<<","<<fakeHit->getV()<<")"<<std::endl;
      delete fakeHit;
      
      // Loop over all hits found and check if any have a sufficiently good delta Chi2
      vector<KDCluster*> goodHits;
      KDCluster* bestCluster = NULL;
      
//      std::cout<<"- track "<<currentTrack<<" has "<<results2.size()<<" nearest neighbours"<<std::endl;
//      for(int newHit=0;newHit<results2.size();newHit++){
      unsigned int nKDHits = tempClusterContainer.size();
      for(unsigned int nKDHit = 0; nKDHit<nKDHits; nKDHit++){
        
        // Get the kdHit and check if it has already been used (assigned to a track)
        KDCluster* kdhit = tempClusterContainer[nKDHit];
        if( kdhit->used() )continue;
        
        // First check that the hit is not wildly away from the track (make cell and check angle)
//        Cell* extensionCell = new Cell(conformalTracks[currentTrack]->m_clusters[0],results2[newHit]);
        Cell* extensionCell = new Cell(conformalTracks[currentTrack]->m_clusters[0],kdhit);
        double cellAngle = seedCell->getAngle(extensionCell);
        double cellAngleRZ = seedCell->getAngleRZ(extensionCell);
        delete extensionCell;
        if( cellAngle > 3.*m_maxCellAngle || cellAngleRZ > 3.*m_maxCellAngleRZ ){
          continue;
        }

        // Now fit the track with the new hit and check the increase in chi2
//        double deltaChi2 = fitWithPoint(*conformalTracks[currentTrack],results2[newHit]); //conformalTracks[currentTrack]->deltaChi2(results2[newHit]);
//        double deltaChi2 = fitWithPoint(*(conformalTracks[currentTrack]->kalmanTrack()),kdhit); //conformalTracks[currentTrack]->deltaChi2(results2[newHit]);
        double deltaChi2 = fitWithPoint(*(conformalTracks[currentTrack]),kdhit); //conformalTracks[currentTrack]->deltaChi2(results2[newHit]);
//        std::cout<<"- delta chi2 of hit "<<newHit<<" is "<<deltaChi2<<". Hit position ("<<results2[newHit]->getU()<<","<<results2[newHit]->getV()<<")"<<std::endl;
        if(deltaChi2 > 100.) continue;
        std::cout<<"Good hit on track " <<currentTrack<<std::endl;
//        results2[newHit]->setDeltaChi2(deltaChi2);
//        goodHits.push_back(results2[newHit]);
        
        // Check if there is a best cluster on this layer
        if(bestCluster == NULL){
          std::cout<<"Replacing null best cluster"<<std::endl;
          bestCluster = kdhit;
          bestCluster->setDeltaChi2(deltaChi2);
        }else{
          // If new layer, save the old best cluster
          if(!kdhit->sameLayer(bestCluster)){
            goodHits.push_back(bestCluster);
            bestCluster->used(true);
            
            // Recalculate chi2 with new hit on kalman track
            
            
            kdhit->setDeltaChi2(deltaChi2);
            bestCluster = kdhit;
            // If same layer, replace best cluster if chi2 is better
          }else if(deltaChi2 < bestCluster->getDeltaChi2()){
            kdhit->setDeltaChi2(deltaChi2);
            bestCluster = kdhit;
            std::cout<<"Replacing best cluster with higher chi2"<<std::endl;
          }
        }
        
        
        // Now we have a point to consider on the current layer. If we are moving onto a new layer then take the best hit and attach it
//        if(nKDHit == (nKDHits-1) || !(kdhit->sameLayer(tempClusterContainer[nKDHit+1]))){
//          
//          std::cout<<"Last cluster/next hit on different layer"<<std::endl;
//
//          if(bestCluster == NULL) continue;
//          
//          std::cout<<"Cluster push back"<<std::endl;
//
////          conformalTracks[currentTrack]->kalmanTrack()->addCluster(bestCluster);
//          goodHits.push_back(bestCluster);
//          bestCluster->used(true);
//          bestCluster = NULL;
//
//        }

//        kdhit->setDeltaChi2(deltaChi2);
//        goodHits.push_back(kdhit);
//        conformalTracks[currentTrack]->kalmanTrack()->addCluster(kdhit);
//        kdhit->used(true);
      }
      
      if(bestCluster != NULL){
        std::cout<<"Cluster post-push back"<<std::endl;
        goodHits.push_back(bestCluster);
        bestCluster->used(true);
        bestCluster = NULL;
      }
      
      delete seedCell;
      
      for(int i=0;i<goodHits.size();i++){
        conformalTracks[currentTrack]->add(goodHits[i]);
        goodHits[i]->used(true);
      }
//      if(goodHits.size() > 0){
//        conformalTracks[currentTrack]->linearRegression();
//        conformalTracks[currentTrack]->linearRegressionConformal();
//      }
      std::cout<<"- pushed back "<<goodHits.size()<<" good hits to track "<<currentTrack<<std::endl;
      
    }
 
  	// ---------------------------------------------------------------------
		// Try to create new tracks using all of the kdHits currently held
	  // ---------------------------------------------------------------------
		unsigned int nKDHits = kdClusters.size();
    streamlog_out( DEBUG4 )<<"Seeding with hits"<<std::endl;
    streamlog_out( DEBUG4 )<<"Attempting to seed with hits: "<<nKDHits<<std::endl;
    
    // Loop over all current hits, using only vertex detectors as seeds
    if(collection > 1) continue;
//    std::cout<<"Seeding with hits. Max collection number "<<collection<<", containing a total of "<<nKDHits<<" hits"<<std::endl;
    for(unsigned int nKDHit = 0; nKDHit<nKDHits; nKDHit++){
      
      // Get the kdHit and check if it has already been used (assigned to a track)
      KDCluster* kdhit = kdClusters[nKDHit];
//      std::cout<<"Seeding with hit "<<nKDHit<<std::endl;
      if(debugSeed && kdhit == debugSeed) std::cout<<"Starting to seed with debug cluster"<<std::endl;
//      if( kdhit->used() ) continue;
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
      double theta = kdhit->getTheta();
      nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);
//      nearestNeighbours->allNeighboursInRadius(kdhit, m_maxDistance, results);
      
      // Sort the neighbours from outer to inner radius
      if(debugSeed && kdhit == debugSeed) std::cout<<"- picked up "<<results.size()<<" neighbours from theta search"<<std::endl;
      if(results.size() == 0) continue;
      std::sort(results.begin(),results.end(),sort_by_radiusKD);
      
      // Objects to hold cells
      std::vector<Cell*> cells;
      
      // Make seed cells pointing inwards (conformal space)
      for(unsigned int neighbour=0;neighbour<results.size();neighbour++){
        
        // Get the neighbouring hit
        KDCluster* nhit = results[neighbour];
        
        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
//        if( nhit->used() ) continue;
        if(kdhit->sameLayer(nhit)) continue;
        if(nhit->getR() >= kdhit->getR()) continue;
        
        // Check if the cell would be too long (hit very far away)
        double length = sqrt( (kdhit->getU()-nhit->getU())*(kdhit->getU()-nhit->getU()) + (kdhit->getV()-nhit->getV())*(kdhit->getV()-nhit->getV()) );
        if(length > m_maxDistance) continue;

        // Create the new seed cell
        Cell* cell = new Cell(kdhit,nhit);
        cells.push_back(cell);
        if(debugSeed && kdhit == debugSeed) std::cout<<"- made cell with neighbour "<<neighbour<<" at "<<nhit->getU()<<","<<nhit->getV()<<std::endl;
        
        // Debug plotting
        if(m_debugPlots && m_eventNumber == 0 && collection == 2){
          m_canvConformalEventDisplayAllCells->cd();
          drawline(kdhit,nhit,1);
        }

      }
      
      if(debugSeed && kdhit == debugSeed) std::cout<<"- produced "<<cells.size()<<" seed cells"<<std::endl;
      
      // No seed cells produced
      if(cells.size() == 0) continue;
      
      // All seed cells have been created, now try create all "downstream" cells until no more can be added
      extendSeedCells(cells, used, nearestNeighbours, false, debugHits);

      // Now have all cells stemming from this seed hit. If it is possible to produce a track (ie. cells with depth X) then we will now...
      //      if(depth < (m_minClustersOnTrack-1)) continue; // TODO: check if this is correct
      
      // We create all acceptable tracks by looping over all cells with high enough weight to create
      // a track and trace their route back to the seed hit. We then have to choose the best candidate
      // at the end (by minimum chi2 of a linear fit)
      std::map<Cell*,bool> usedCells;
      std::map<Cell*,bool> usedCells2;
      std::vector<KDTrack*> cellTracks;
      
      // Sort Cells from highest to lowest weight
      std::sort(cells.begin(),cells.end(),sort_by_cellWeight);

      // Create tracks by following a path along cells
      int nCells = cells.size();
      for(int itCell=0;itCell<nCells;itCell++){
        
        // Check if this cell has already been used
        if(debugSeed && kdhit == debugSeed) std::cout<<"-- looking at cell "<<itCell<<std::endl;
        if(usedCells.count(cells[itCell])) continue;

        // Check if this cell could produce a track (is on a long enough chain)
        if(cells[itCell]->getWeight() < (m_minClustersOnTrack-2)) break;
        
        // Produce all tracks leading back to the seed hit from this cell
        std::vector<cellularTrack*> candidateTracks;
        std::vector<cellularTrack*> candidateTracksTemp;
        createTracksNew(candidateTracksTemp,cells[itCell],usedCells2); // Move back to using used cells here? With low chi2/ndof?
        
        for(int itTrack=0;itTrack<candidateTracksTemp.size();itTrack++){
          if(candidateTracksTemp[itTrack]->size() >= (m_minClustersOnTrack-1) ) candidateTracks.push_back(candidateTracksTemp[itTrack]);
          else delete candidateTracksTemp[itTrack];
        }

        // Debug plotting
        if(m_debugPlots && m_eventNumber == 2){
          m_canvConformalEventDisplayAcceptedCells->cd();
          for(int iTr=0;iTr<candidateTracks.size();iTr++){
            cellularTrack* track = candidateTracks[iTr];
            for(unsigned int trackCell=0;trackCell<track->size();trackCell++){
              drawline((*track)[trackCell]->getStart(),(*track)[trackCell]->getEnd(),track->size()-trackCell);
            }
          }
        }

        // Look at the candidate tracks and fit them (+pick the best chi2/ndof)
//        std::cout<<"- produced "<<candidateTracks.size()<<" candidate tracks"<<std::endl;
        if(debugSeed && kdhit == debugSeed) std::cout<<"- produced "<<candidateTracks.size()<<" candidate tracks"<<std::endl;
        if(candidateTracks.size() == 0) continue;
        std::vector<double> chi2ndof;
        std::vector<KDTrack*> bestTracks;

        getFittedTracks(bestTracks,candidateTracks,usedCells); // Returns all tracks at the moment, not lowest chi2 CHANGE ME

        // Store track(s) for later
        cellTracks.insert(cellTracks.end(),bestTracks.begin(),bestTracks.end());
        
      }
      
      // All tracks leading back to the seed hit have now been found. Decide which are the feasible candidates (may be more than 1)
      if(debugSeed && kdhit == debugSeed) std::cout<<"== final number of candidate tracks to this seed hit: "<<cellTracks.size()<<std::endl;
      if(cellTracks.size() == 0){
        // Clean up
        for(unsigned int itCell=0;itCell<cells.size();itCell++) delete cells[itCell];
        continue;
      }

      std::vector<KDTrack*> bestTracks = cellTracks; //CHANGE ME - temp to give all tracks
//      std::vector<KDTrack*> bestTracks = getLowestChi2(cellTracks,cellTracksChi2ndof,chi2ndof);
//      std::cout<<"== final number of stored tracks to this seed hit: "<<bestTracks.size()<<std::endl;
      if(debugSeed && kdhit == debugSeed){
        std::cout<<"== final number of stored tracks to this seed hit: "<<bestTracks.size()<<std::endl;
        for(int itBest=0;itBest<bestTracks.size();itBest++) std::cout<<"- track "<<itBest<<" has chi2/ndof "<<bestTracks[itBest]->chi2ndof()<<std::endl;
      }

      // Could now think to do the full helix fit and apply a chi2 cut. TODO
  
      // Store the final CA tracks. First sort them by length, so that if we produce a clone with just 1 hit missing, the larger track is taken
      std::sort(bestTracks.begin(),bestTracks.end(),sort_by_length);
      for(unsigned int itTrack=0;itTrack<bestTracks.size();itTrack++){
        bool bestTrackUsed = false;
       
        // Cut on chi2
//        if(collection != trackerHitCollections.size() && (bestTracks[itTrack]->chi2ndof() > m_chi2cut || bestTracks[itTrack]->chi2ndofZS() > m_chi2cut)) {
        if(bestTracks[itTrack]->chi2ndof() > m_chi2cut || bestTracks[itTrack]->chi2ndofZS() > m_chi2cut) {
          delete bestTracks[itTrack];
          bestTracks[itTrack] = NULL;
          continue;
        }
        
        // Check if the new track is a clone
        bool clone = false;
      
        for(int existingTrack=0;existingTrack<conformalTracks.size();existingTrack++){
          
          const int nOverlappingHits = overlappingHits(bestTracks[itTrack],conformalTracks[existingTrack]);
          if( nOverlappingHits >= 2) {
            clone = true;
            
            // Calculate the new and existing chi2 values
            double newchi2 = (bestTracks[itTrack]->chi2ndofZS()*bestTracks[itTrack]->chi2ndofZS() + bestTracks[itTrack]->chi2ndof()*bestTracks[itTrack]->chi2ndof());
            double oldchi2 = (conformalTracks[existingTrack]->chi2ndofZS()*conformalTracks[existingTrack]->chi2ndofZS() + conformalTracks[existingTrack]->chi2ndof()*conformalTracks[existingTrack]->chi2ndof());
            
            double deltachi2ZS = (bestTracks[itTrack]->chi2ndofZS()-conformalTracks[existingTrack]->chi2ndofZS());
            double deltachi2 = (bestTracks[itTrack]->chi2ndof()-conformalTracks[existingTrack]->chi2ndof());
            
            // If the new track is an existing track + segment, make sure the increase in chi2 is not too much
            if(nOverlappingHits == conformalTracks[existingTrack]->m_clusters.size()){
              
              // Increase in chi2 is too much (double)
              if( (newchi2 - oldchi2) > oldchi2) break;
              
              // Otherwise replace the existing track
              delete conformalTracks[existingTrack];
              conformalTracks[existingTrack] = bestTracks[itTrack];
              bestTrackUsed = true;
            }
            
            // If the new track is a subtrack of an existing track, don't consider it further (already try removing bad hits from tracks
            else if(nOverlappingHits == bestTracks[itTrack]->m_clusters.size()) break;
            
            // Otherwise take the longest if the delta chi2 is not too much
            else if(bestTracks[itTrack]->m_clusters.size() >= conformalTracks[existingTrack]->m_clusters.size()){ // New track longer/equal in length
              
              // Increase in chi2 is too much (double)
              if( (newchi2 - oldchi2) > oldchi2) break;
              
              // Otherwise take it
              delete conformalTracks[existingTrack];
              conformalTracks[existingTrack] = bestTracks[itTrack];
              bestTrackUsed = true;
            }
            else if(bestTracks[itTrack]->m_clusters.size() < conformalTracks[existingTrack]->m_clusters.size()){ // Old track longer
              
              // Must improve chi2 by factor two
              if( (newchi2 - 0.5*oldchi2) > 0.) break;
              
              // Otherwise take it
              delete conformalTracks[existingTrack];
              conformalTracks[existingTrack] = bestTracks[itTrack];
              bestTrackUsed = true;
            }
            
            break;
          }
        }
      
        // If not a clone, save the new track
        if(!clone){
          conformalTracks.push_back(bestTracks[itTrack]);
          bestTrackUsed = true;
          if(debugSeed && kdhit == debugSeed)
            std::cout<<"== Pushing back best track with chi2/ndof "<<bestTracks[itTrack]->chi2ndof()<<std::endl;
          else
            std::cout<<"Pushing back best track with chi2/ndof "<<bestTracks[itTrack]->chi2ndof()<<std::endl;
        }
        
        if( not bestTrackUsed ) {
          delete bestTracks[itTrack];
          bestTracks[itTrack] = NULL;
        }
        
      } // end for besttracks
      
      // Clean up
      for(unsigned int itCell=0;itCell<cells.size();itCell++) delete cells[itCell];
      
    }
    // Now finished looking at this hit collection. All tracks which can have extra hits added
    // have had them added, and no new tracks were possible using the sum of all collections
    // till now. Add the next collection and try again...
    delete nearestNeighbours;
    
    // Mark hits from "good" tracks as being used
    for(unsigned int itTrack=0;itTrack<conformalTracks.size();itTrack++){
        for(unsigned int itHit=0;itHit<conformalTracks[itTrack]->m_clusters.size();itHit++){
          conformalTracks[itTrack]->m_clusters[itHit]->used(true);
        }
    }

  }
  
    
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
  std::cout<<"*** CA has made "<<conformalTracksFinal.size()<< (conformalTracksFinal.size() == 1 ? " track ***" : " tracks ***") <<std::endl;

  // Loop over all track candidates
  for(unsigned int caTrack=0;caTrack<conformalTracksFinal.size();caTrack++){
    
    // Vector of all the hits on the track
    KDTrack* conformalTrack = conformalTracksFinal[caTrack];
    streamlog_out( DEBUG5 )<<"Made a track with "<<conformalTrack->m_clusters.size()<<" hits"<<std::endl;
    
    // Make the LCIO track hit vector
    EVENT::TrackerHitVec trackHits;
    for(unsigned int itHit=0;itHit<conformalTrack->m_clusters.size();itHit++){
      KDCluster* cluster = conformalTrack->m_clusters[itHit];
      trackHits.push_back(kdClusterMap[cluster]);
    }
    // Add kalman filtered hits
    if(conformalTrack->kalmanTrack() != NULL){
      
      KalmanTrack* kalmanTrack = conformalTrack->kalmanTrack();
      for(int i=1;i<kalmanTrack->m_kalmanClusters.size();i++){
        KDCluster* cluster = kalmanTrack->m_kalmanClusters[i];
        trackHits.push_back(kdClusterMap[cluster]);
      }
    }
    
    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(),trackHits.end(),sort_by_radius);
    
    // Now we can make the track object and relations object, and fit the track
    TrackImpl* track = new TrackImpl ;
    
    // First, for some reason there are 2 track objects, one which gets saved and one which is used for fitting. Don't ask...
    MarlinTrk::IMarlinTrack* marlinTrack = trackFactory->createTrack();
    
    // Make an initial covariance matrix with very broad default values
    EVENT::FloatVec covMatrix (15,0); // Size 15, filled with 0s
    covMatrix[0]  = ( m_initialTrackError_d0    ); //sigma_d0^2
    covMatrix[2]  = ( m_initialTrackError_phi0  ); //sigma_phi0^2
    covMatrix[5]  = ( m_initialTrackError_omega ); //sigma_omega^2
    covMatrix[9]  = ( m_initialTrackError_z0    ); //sigma_z0^2
    covMatrix[14] = ( m_initialTrackError_tanL  ); //sigma_tanl^2
    
    // Try to fit
//    int fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackHits, track, MarlinTrk::IMarlinTrack::forward, covMatrix, m_magneticField, m_maxChi2perHit);
    
    // Check track quality - if fit fails chi2 will be 0. For the moment add hits by hand to any track that fails the track fit, and store it as if it were ok...
//    if(track->getChi2() == 0.){
//      std::cout<<"Fit failed. Track has "<<track->getTrackerHits().size()<<" hits"<<std::endl;
//      std::cout<<"Fit fail error "<<fitError<<std::endl;
      for(unsigned int p=0;p<trackHits.size();p++){
        track->addHit(trackHits[p]);
      }
//    }//
    delete marlinTrack;
//    delete track; delete marlinTrack; continue;}
    
    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ] = trackHits.size();
    
    // Push back to the output container
    trackCollection->addElement(track);
    
  }
  
  // calculate purities and check if tracks have been reconstructed
  if(m_debugPlots){

    for(int itrack=0; itrack<conformalTracksFinal.size();itrack++){
      KDTrack* debugTrack = conformalTracksFinal[itrack];
      
      m_conformalChi2->Fill(debugTrack->chi2ndof());
      double purity = checkReal(debugTrack,kdParticles,reconstructed,particleHits);
      if(purity >= m_purity){
        m_conformalChi2real->Fill(debugTrack->chi2ndof());
      }
      if(purity < m_purity){
        m_conformalChi2fake->Fill(debugTrack->chi2ndof());
      }
      m_conformalChi2Purity->Fill(purity,debugTrack->chi2ndof());
    }
  }//*/
  
  // Draw the cells for all produced tracks
  if(m_debugPlots && m_eventNumber == 0){
    m_canvConformalEventDisplay->cd();
    for(int itrack=0; itrack<conformalTracksFinal.size();itrack++){
      KDTrack* debugTrack = conformalTracksFinal[itrack];
      std::vector<KDCluster*> clusters = debugTrack->m_clusters;
      for(int itCluster=1;itCluster<clusters.size();itCluster++) drawline(clusters[itCluster-1],clusters[itCluster],clusters.size()-itCluster);
    }
  }
  
  // Draw the conformal event display hits for debugging
  if(m_debugPlots && m_eventNumber == 0){
    m_canvConformalEventDisplay->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayAllCells->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayAcceptedCells->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayMCunreconstructed->cd();
    m_conformalEvents->DrawCopy("same");
  }
  if(m_debugPlots){
    int nReconstructed(0), nUnreconstructed(0);
    // Additionally draw all tracks that were not reconstructed
    m_canvConformalEventDisplayMCunreconstructed->cd();
    int nParticles = particleCollection->getNumberOfElements();
    
    // Make the nearest neighbour tree to debug particle reconstruction issues
    vector<KDCluster*> kdClusters;
    for(int i=0;i<trackerHitCollections.size();i++) kdClusters.insert(kdClusters.begin(),collectionClusters[i].begin(),collectionClusters[i].end());
    KDTree* nearestNeighbours = new KDTree(kdClusters,m_thetaRange);
    
    for(int itP=0;itP<nParticles;itP++){
      // Get the particle
      MCParticle* mcParticle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itP) ) ;
      // Get the conformal hits
      if(particleHits.count(mcParticle) == 0) continue;
      std::vector<KDCluster*> mcHits = particleHits[mcParticle];
      // Cut on the number of hits
      int uniqueHits = getUniqueHits(mcHits);
      if(uniqueHits < m_minClustersOnTrack) continue;
      // Check if it was stable
      if(mcParticle->getGeneratorStatus() != 1) continue;
      // Check if it was reconstructed
      if(reconstructed.count(mcParticle)){nReconstructed++; continue;}
      // Draw the cells connecting the hits
      std::sort(mcHits.begin(),mcHits.end(),sort_by_radiusKD);
      for(int iHit=0;iHit<(mcHits.size()-1);iHit++){
        drawline(mcHits[iHit],mcHits[iHit+1],iHit+1);
      }
      // List the pt of unreconstructed particles
      double particlePt = sqrt( mcParticle->getMomentum()[0]*mcParticle->getMomentum()[0] + mcParticle->getMomentum()[1]*mcParticle->getMomentum()[1] );
      std::cout<<"Unreconstructed particle pt: "<<particlePt<<std::endl;
      nUnreconstructed++;
      
      // Check why particles were not reconstructed
      checkReconstructionFailure(mcParticle, particleHits, used, nearestNeighbours);
    }
    std::cout<<"Reconstructed "<<nReconstructed<<" particles out of "<<nReconstructed+nUnreconstructed<<". Gives efficiency "<<100.*(double)nReconstructed/(double)(nReconstructed + nUnreconstructed)<<"%"<<std::endl;
    delete nearestNeighbours;
  }
  
  // Save the output track collection
  evt->addCollection( trackCollection , m_outputTrackCollection ) ;
  evt->addCollection( debugHitCollection, m_outputDebugHits );
  
  // Increment the event number
  m_eventNumber++ ;
  
}

void ConformalTracking::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void ConformalTracking::end(){
  
  streamlog_out(MESSAGE) << " end()  " << name()
  << " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
  << std::endl ;
  
  // Write debug canvases to output file
  if(m_debugPlots){
    m_canvConformalEventDisplay->Write();
    m_canvConformalEventDisplayAllCells->Write();
    m_canvConformalEventDisplayAcceptedCells->Write();
    m_canvConformalEventDisplayMC->Write();
    m_canvConformalEventDisplayMCunreconstructed->Write();
  }
  
  //FIXME trackFactory is leaking Memory, but probably a MarlinTRK issue

}

// Get a collection from the event object
void ConformalTracking::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG4 )<< "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

// Extend seed cells
void ConformalTracking::extendSeedCells(std::vector<Cell*>& cells, std::map<KDCluster*,bool> used, KDTree* nearestNeighbours, bool extendingTrack, std::vector<KDCluster*> debugHits){
  
  unsigned int nCells=0; int depth = 0; int startPos=0;
  
  // Keep track of existing cells in case there are branches in the track
  std::map<KDCluster*, std::vector<Cell*> > existingCells;
  
  // Try to create all "downstream" cells until no more can be added
  while(cells.size() != nCells){
    
    // Extend all cells with depth N. In the next iteration, look at cells with depth N+1
    nCells = cells.size();
    for(unsigned int itCell=startPos;itCell<nCells;itCell++){
      
      // Get the end point of the cell (to search for neighbouring hits to form new cells connected to this one)
      KDCluster* hit = cells[itCell]->getEnd();
      double searchDistance = m_maxDistance; //hit->getR();
//      if(searchDistance > hit->getR()) searchDistance = 1.2*hit->getR();
      
      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      KDCluster* fakeHit = extrapolateCell(cells[itCell],searchDistance/2.); // TODO: make this search a function of radius
      VecCluster results;
      nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25*searchDistance/2., results);
      delete fakeHit;
      
      if(extendingTrack) std::cout<<"Found "<<results.size()<<" neighbours from cell extrapolation"<<std::endl;
      // Make new cells pointing inwards
      for(unsigned int neighbour=0;neighbour<results.size();neighbour++){
        
        if(extendingTrack) std::cout<<"looking at neighbour "<<neighbour<<std::endl;
        // Get the neighbouring hit
        KDCluster* nhit = results[neighbour];
        
        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
//        if(used.count(nhit)){if(extendingTrack) std::cout<<"- used"<<std::endl; continue;}
//        if( nhit->used() )continue;
        if(hit->sameLayer(nhit)){if(extendingTrack) std::cout<<"- same layer"<<std::endl; continue;}
        if(nhit->getR() >= hit->getR()){if(extendingTrack) std::cout<<"- higher radius"<<std::endl; continue;}
//        double zdifference = abs(hit->getZ()-nhit->getZ());
//        if(zdifference > 10.){
//          if(hit->getZ() > 0. && nhit->getZ() < hit->getZ()){if(extendingTrack) std::cout<<"- z cut"<<std::endl; continue;}
//          if(hit->getZ() < 0. && nhit->getZ() > hit->getZ()){if(extendingTrack) std::cout<<"- z cut"<<std::endl; continue;}
//        }
        
        // Check if this cell already exists (rejoining branch) FIXME - allows rejoining a branch without checking cell angles
        if(existingCells.count(hit) != 0){
          bool alreadyExists=false;
          int nExistingCells = existingCells[hit].size();
          for(int iterExisting=0;iterExisting<nExistingCells;iterExisting++){
            if( existingCells[hit][iterExisting]->getEnd() == nhit ){
              alreadyExists = true;
              
              // Check if cell angle is too large to rejoin
              if( cells[itCell]->getAngle(existingCells[hit][iterExisting]) > m_maxCellAngle || cells[itCell]->getAngleRZ(existingCells[hit][iterExisting]) > m_maxCellAngleRZ ) continue;
              
              // Otherwise add the path
              cells[itCell]->setTo(existingCells[hit][iterExisting]);
              existingCells[hit][iterExisting]->setFrom(cells[itCell]);
              updateCell(existingCells[hit][iterExisting]);
            }
          }
          if(alreadyExists) continue;
        }
        
        // Make the new cell
        Cell* cell = new Cell(hit,nhit);
        
        // Check if the new cell is compatible with the previous cell (angle between the two is acceptable)
//        if( cells[itCell]->getAngle(cell) > (m_maxCellAngle*exp(-0.001/nhit->getR())) ){
        if( cells[itCell]->getAngle(cell) > m_maxCellAngle || cells[itCell]->getAngleRZ(cell) > m_maxCellAngleRZ ){
          
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
    startPos=nCells;
    depth++;
  }
  
  // No more downstream cells can be added
  
}

// Draw a line on the current canvas
void ConformalTracking::drawline(KDCluster* hitStart, KDCluster* hitEnd, int colour, int style){
  
  TLine *line = new TLine(hitStart->getU(),hitStart->getV(),hitEnd->getU(),hitEnd->getV());
  line->SetLineColor(colour);
  line->SetLineStyle(style);
  line->Draw();
  
}

// New test at creating cellular tracks. In this variant, don't worry about clones etc, give all possible routes back to the seed cell. Then cut
// on number of clusters on each track, and pass back (good tracks to then be decided based on best chi2
void ConformalTracking::createTracksNew(std::vector<cellularTrack*>& finalcellularTracks, Cell* seedCell, std::map<Cell*,bool>& usedCells){
 
 	// Final container to be returned
  std::vector<cellularTrack*> cellularTracks;
  
  // Make the first cellular track using the seed cell
  cellularTrack* seedTrack = new cellularTrack();
  cellularTracks.push_back( seedTrack );
  seedTrack->push_back(seedCell);
  
  // Now start to follow all paths back from this seed cell
  // While there are still tracks that are not finished (last cell weight 0), keep following their path
  while(toBeUpdated(cellularTracks)){
    
//    std::cout<<"== Updating "<<cellularTracks.size()<<" tracks"<<std::endl;
    // Loop over all (currently existing) tracks
    int nTracks = cellularTracks.size();
    for(int itTrack=0;itTrack<nTracks;itTrack++){
     
      // If the track is finished, do nothing
//      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if(cellularTracks[itTrack]->back()->getFrom()->size() == 0){
//        std::cout<<"-- Track "<<itTrack<<" is finished"<<std::endl;
        continue;
      }
      
      // While there is only one path leading from this cell, follow that path
      Cell* cell = cellularTracks[itTrack]->back();
//      std::cout<<"-- Track "<<itTrack<<" has "<<(*(cell->getFrom())).size()<<" cells attached to the end of it"<<std::endl;
//      while(cell->getWeight() > 0 && (*(cell->getFrom())).size() == 1){
      while((*(cell->getFrom())).size() == 1){
//        std::cout<<"- simple extension"<<std::endl;
        // Get the cell that it attaches to
        Cell* parentCell = (*(cell->getFrom()))[0];
        // Attach it to the track and continue
        cellularTracks[itTrack]->push_back(parentCell);
        cell = parentCell;
      }

      // If the track is finished, do nothing
//      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if(cellularTracks[itTrack]->back()->getFrom()->size() == 0) continue;

      // If the weight is != 0 and there is more than one path to follow, branch the track (create a new one for each path)
      int nBranches = (*(cell->getFrom())).size();
      
//      std::cout<<"- making "<<nBranches<<" branches"<<std::endl;

      // For each additional branch make a new track
      for(int itBranch=1;itBranch<nBranches;itBranch++){
        cellularTrack* branchedTrack = new cellularTrack();
        (*branchedTrack) = (*cellularTracks[itTrack]);
        cellularTracks.push_back( branchedTrack );
        branchedTrack->push_back((*(cell->getFrom()))[itBranch]);
      }
      
      // Keep the existing track for the first branch
      cellularTracks[itTrack]->push_back((*(cell->getFrom()))[0]);
      
    }
    
  }

  int nTracks = cellularTracks.size();
  for(int itTrack=0;itTrack<nTracks;itTrack++){
//    if(cellularTracks[itTrack]->size() >= (m_minClustersOnTrack-1) )
      finalcellularTracks.push_back(cellularTracks[itTrack]);
  }
  
  return;
  
}

// Check if any of the tracks in a collection still have to be updated
bool ConformalTracking::toBeUpdated(std::vector<cellularTrack*>const& cellularTracks){
  bool update=false;
  for(unsigned int iTrack=0;iTrack<cellularTracks.size();iTrack++) if( cellularTracks[iTrack]->back()->getFrom()->size() > 0 ){update = true; break;}
  return update;
}

// Given a list of connected cells (so-called cellular tracks), return the candidate(s) with lowest chi2/degrees of freedom.
// Attempt to remove hits to see if there is a significant improvement in the chi2/ndof, to protect against noise hits causing
// a good track to be discarded. If several candidates have low chi2/ndof and are not clones (limited sharing of hits) then
// return all of them. Given that there is no material scattering taken into account this helps retain low pt tracks, which
// may have worse chi2/ndof than ghosts/real tracks with an additional unrelated hit from the low pt track.
void ConformalTracking::getFittedTracks(std::vector<KDTrack*>& finalTracks, std::vector<cellularTrack*>& candidateTracks, std::map<Cell*,bool>& usedCells){
  
  // Make a container for all tracks being considered, initialise variables
  std::vector<KDTrack*> trackContainer;
//  std::vector<double> trackChi2ndofs;
  
  // Loop over all candidate tracks and do an inital fit to get the track angle (needed to calculate the
  // hit errors for the error-weighted fit)
  for(unsigned int itTrack=0;itTrack<candidateTracks.size();itTrack++){
    
    // If there are not enough hits on the track, ignore it
    if(candidateTracks[itTrack]->size() < (m_minClustersOnTrack - 2)) {
      delete candidateTracks[itTrack];
      candidateTracks[itTrack] = NULL;
      continue;
    }
    
		// Make the fitting object. TGraphErrors used for 2D error-weighted fitting
    KDTrack* track = new KDTrack();

    // Loop over all hits and add them to the fitter (and track)
    double npoints=0.;
    KDCluster* kdStart = (*candidateTracks[itTrack])[0]->getEnd();
    track->add(kdStart); npoints++;
    
    for(unsigned int trackCell=0;trackCell<(*candidateTracks[itTrack]).size();trackCell++){
      KDCluster* kdEnd = (*candidateTracks[itTrack])[trackCell]->getStart();
      track->add(kdEnd); npoints++;
    }
/*
    // Set up the track fitting
//    std::cout<<"LOOK for me"<<std::endl;
    ROOT::Math::Functor FCNFunction(track,2);
    newFitter.SetFunction(FCNFunction);
//    globalTrack = &track;
    newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
    newFitter.SetVariable(1,"intercept", 0., 0.1);
    
    // Fit the track, first in uv space, then in sz space
    track.setConformalFit(true);
    newFitter.Minimize();
		
    // Now set the track parameters from the conformal fit, and fit in sz
    track.setGradient(newFitter.X()[0]);
    track.setIntercept(newFitter.X()[1]); */
/*    track.setGradientError(newFitter.Errors()[0]);
    track.setInterceptError(newFitter.Errors()[1]);
    track.setConformalFit(false);
    
    double b = 1./(2.*track.intercept());
    double a = -1.*b*track.gradient();
    
    double xMa = track.clusters()[npoints-1]->getX() - a;
    double yMb = track.clusters()[npoints-1]->getY() - b;
    double s = atan2(yMb,xMa);

    double xMa1 = track.clusters()[0]->getX() - a;
    double yMb1 = track.clusters()[0]->getY() - b;
    double s1 = atan2(yMb1,xMa1);

    double startingGuess = (s-s1)/(track.clusters()[npoints-1]->getZ()-track.clusters()[0]->getZ());
    double startingIntercept = s-(startingGuess*track.clusters()[npoints-1]->getZ());
    
    std::cout<<"-- Starting guess for sz gradient is "<<startingGuess<<std::endl;
    std::cout<<"-- Starting guess for sz intercept is "<<startingIntercept<<std::endl;
    newFitter.Clear();
    
    ROOT::Math::Functor FCNFunction2(track,2);
    newFitter.SetFunction(FCNFunction2);

    newFitter.SetVariable(0,"gradient", startingGuess, 0.0001);
    newFitter.SetVariable(1,"intercept", startingIntercept, 0.1);
    newFitter.Minimize();

    track.setGradientZS(newFitter.X()[0]);
    track.setInterceptZS(newFitter.X()[1]);
    //*/
    // Calculate the track chi2 with the final fitted values
//    track.fit();
//    std::cout<<"-- Track fitting gives gradient of "<<newFitter.X()[0]<<", intercept of "<<newFitter.X()[1]<<" and chi2 of "<<track.chi2()<<std::endl;
    track->linearRegression();
    track->linearRegressionConformal(); //FCC study
//    double chi2sz = track.calculateChi2SZ();
//    if(chi2sz > 2.1e-08 && chi2sz < 2.2e-08){
//      track.FillDistribution(m_szDistribution);
//      for(int i=0;i<track.nPoints();i++){
//        m_uvDistribution->Fill(track.clusters()[i]->getU(),track.clusters()[i]->getV());
//        m_xyDistribution->Fill(track.clusters()[i]->getX(),track.clusters()[i]->getY());
//        m_xyzDistribution->Fill(track.clusters()[i]->getX(),track.clusters()[i]->getY(),track.clusters()[i]->getZ());
//      }
//    }
    
//    std::cout<<"Done looking"<<std::endl;
//    double chi2ndof = track->chi2ndofZS();
    double chi2ndof = track->chi2()/(npoints-2); //FCC study
    
//    if(track.calculateChi2SZ() > 1.e6) continue;
    
    // We try to see if there are spurious hits causing the chi2 to be very large. This would cause us to throw away
    // good tracks with perhaps just a single bad hit. Try to remove each hit and see if fitting without it causes a
    // significant improvement in chi2/ndof
    
    // Loop over each hit (starting at the back, since we will use the 'erase' function to get rid of them)
    // and see if removing it improves the chi2/ndof
    double removed=0;
    if(chi2ndof > m_chi2cut && chi2ndof < m_chi2cut){ //CHANGE ME?? Upper limit to get rid of really terrible tracks (temp lower changed from 0 to m_chi2cut)
      for(int point=npoints-1;point>=0;point--){
        
        // Stop if we would remove too many points on the track to meet the minimum hit requirement (or if the track has more
        // than 2 hits removed)
        if((npoints-removed-1) < m_minClustersOnTrack || removed == 2) break;
        
        // Refit the track without this point
        double newChi2ndof = fitWithoutPoint(*track,point);
        
        // If the chi2/ndof is significantly better, remove the point permanently CHANGE ME??
//        if( (chi2ndof - newChi2ndof) > 0 && (chi2ndof - newChi2ndof) > 1. ){
        if( (newChi2ndof - chi2ndof) < chi2ndof ){
          track->remove(point);
          removed++;
          chi2ndof = newChi2ndof;
        }
      }
    }
  
    // Store the track information
    trackContainer.push_back(track);
//    trackChi2ndofs.push_back(chi2ndof);
    
    // LOOK AT ME
//    if(chi2ndof < 10.){
//      for(unsigned int trackCell=0;trackCell<candidateTracks[itTrack]->size();trackCell++) usedCells[(*candidateTracks[itTrack])[trackCell]] = true;
//    }

    delete candidateTracks[itTrack];

  }// end for candidateTracks
  
  // TEMP - return all tracks, don't take lowest chi2
//  finalChi2ndofs = trackChi2ndofs;
//  return trackContainer;
  
  // Now have all sets of conformal tracks and their chi2/ndof. Decide which tracks to send back, ie. the one with
  // lowest chi2/ndof, and possibly others if they are not clones and have similar chi2 value
//  std::vector<KDTrack*> finalTracks;
  getLowestChi2(finalTracks,trackContainer);
  
  // Send back the final set of tracks
  return;
  
}

// Pick the lowest chi2/ndof KDTrack from a list of possible tracks, and additionally return other tracks in the collection with similar chi2/ndof values that don't share many hits
void ConformalTracking::getLowestChi2(std::vector<KDTrack*>& finalTracks, std::vector<KDTrack*> trackContainer){
  
  // Get the lowest chi2/ndof value from the given tracks
//  double lowestChi2ndof = *std::min_element(trackChi2ndofs.begin(),trackChi2ndofs.end());
  KDTrack* lowestChi2ndofTrack = trackContainer[0];
  double lowestChi2ndof = lowestChi2ndofTrack->chi2ndof();
  
  for(unsigned int itTrack=0;itTrack<trackContainer.size();itTrack++){
    if(trackContainer[itTrack]->chi2ndof() < lowestChi2ndof){
      lowestChi2ndof = trackContainer[itTrack]->chi2ndof();
      lowestChi2ndofTrack = trackContainer[itTrack];
    }
  }
  
  
	// Final track storage
//  std::vector<KDTrack*> finalTracks;
  
  // Loop over all other tracks and decide whether or not to save them
  for(unsigned int itTrack=0;itTrack<trackContainer.size();itTrack++){
    
    // Look at the difference in chi2/ndof - we want to keep tracks with similar chi2/ndof. If they
    // are clones then take the longest
    if( (trackContainer[itTrack]->chi2ndof() - lowestChi2ndof) < 10. ){
      
      // If same track and longer
//      if(sameTrack(trackContainer[itTrack], lowestChi2ndofTrack)){
//        if(trackContainer[itTrack].nPoints() > lowestChi2ndofTrack.nPoints()){
//        	lowestChi2ndofTrack = trackContainer[itTrack];
//        	lowestChi2ndof = trackChi2ndofs[itTrack];
//        }
//        continue;
//      }
      
      // Store this track
      finalTracks.push_back(trackContainer[itTrack]);

    } else {
      delete trackContainer[itTrack];
      trackContainer[itTrack] = NULL;
    }
  }
  
  // Save the track with the lowest chi2/ndof
//	finalTracks.insert(finalTracks.begin(),lowestChi2ndofTrack);
//	finalChi2ndofs.insert(finalChi2ndofs.begin(),lowestChi2ndof);

  return;
  
}

// Function to check if two KDtracks contain several hits in common
int ConformalTracking::overlappingHits(const KDTrack* track1, const KDTrack* track2){
  
  // Loop over all hits on track 1 and check if that hit is in track 2
  int nHitsInCommon = 0;
//  std::vector<KDCluster*> track1hits = track1->m_clusters;
//  std::vector<KDCluster*> track2hits = track2->m_clusters;
  
  for(int hit=0;hit<track1->m_clusters.size();hit++){
    if( std::find(track2->m_clusters.begin(),track2->m_clusters.end(),track1->m_clusters[hit]) !=  track2->m_clusters.end()) nHitsInCommon++;
  }
  
  // Cut on number of shared hits
//  if(nHitsInCommon > 0.4 * track1hits.size()) return true;
  return nHitsInCommon;
  
}

double ConformalTracking::fitWithoutPoint(KDTrack track,int point){

  // Remove the given point from the track
  track.remove(point);
/*
  // Set up the fitter
  int npoints = track.nPoints();
      ROOT::Math::Functor FCNFunction(track,2);
      newFitter.SetFunction(FCNFunction);
//  globalTrack = &track;
  newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
  newFitter.SetVariable(1,"intercept", 0., 0.1);
  
  // Fit the track, first in uv space, then in sz space
  track.setConformalFit(true);
  newFitter.Minimize();
		
  // Now set the track parameters from the conformal fit, and fit in sz
  track.setGradient(newFitter.X()[0]);
  track.setIntercept(newFitter.X()[1]); */
//  track.setGradientError(newFitter.Errors()[0]);
//  track.setInterceptError(newFitter.Errors()[1]);
//  track.setConformalFit(false);
//  newFitter.Minimize();
//  track.setGradientZS(newFitter.X()[0]);
//  track.setInterceptZS(newFitter.X()[1]);
  
  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal(); // FCC study
//  track.fit();
  double chi2ndof = track.chi2ndof();
  double chi2ndofZS = track.chi2ndofZS(); // FCC study
  
  return sqrt(chi2ndof*chi2ndof + chi2ndofZS*chi2ndofZS) ;
}

void ConformalTracking::updateCell(Cell* cell){
  
  if((*(cell->getTo())).size() != 0){
    for(unsigned int i=0;i<(*(cell->getTo())).size();i++){
      (*(cell->getTo()))[i]->update(cell);
      updateCell((*(cell->getTo()))[i]);
    }
  }
  
}

// Function to extrapolate along a cell in conformal space, producing a fake hit
// a given distance away from the cell endpoint
KDCluster* ConformalTracking::extrapolateCell(Cell* cell, double distance){
  
  // Fake cluster to be returned
  KDCluster* extrapolatedCluster = new KDCluster();
  
  // Linear extrapolation of the cell - TODO: check that cell gradients have correct sign and remove checks here
  double gradient = cell->getGradient();
  double deltaU = sqrt( distance*distance / (1+gradient*gradient) );
  double deltaV = abs(gradient)*deltaU;
  
  if( (cell->getStart()->getU() - cell->getEnd()->getU()) > 0 ) deltaU *= (-1.);
  if( (cell->getStart()->getV() - cell->getEnd()->getV()) > 0 ) deltaV *= (-1.);
  
  extrapolatedCluster->setU( cell->getEnd()->getU() + deltaU );
  extrapolatedCluster->setV( cell->getEnd()->getV() + deltaV );
  
  return extrapolatedCluster;
}
/*
void ConformalTracking::extendTrack(KDTrack* track,std::vector<cellularTrack*> trackSegments, std::map<KDCluster*,bool>& used, std::map<Cell*,bool>& usedCells){

  cout<<"== extending track, have "<<trackSegments.size()<<" candidates"<<endl;
  // For each track segment, perform a kalman filter on the addition points and chose the track extension with the
  // best delta chi2.
  KalmanTrack* bestTrack = NULL; double bestChi2=0.; double bestNpoints;
  for(int nTrackExtension=0;nTrackExtension<trackSegments.size();nTrackExtension++){

    double npoints=0;
//    std::cout<<"-> extension "<<nTrackExtension<<" has "<<trackSegments[nTrackExtension]->size()<<" cells"<<std::endl;
//    for(int c=0;c<trackSegments[nTrackExtension]->size();c++){
//      std::cout<<"-> cell "<<c<<" start ("<<(*trackSegments[nTrackExtension])[c]->getStart()->getU()<<","<<(*trackSegments[nTrackExtension])[c]->getStart()->getV()<<") and stop ("<<(*trackSegments[nTrackExtension])[c]->getEnd()->getU()<<","<<(*trackSegments[nTrackExtension])[c]->getEnd()->getV()<<")"<<std::endl;
//    }
    if( trackSegments[nTrackExtension]->size() < 2) {
      delete trackSegments[nTrackExtension];
      continue;
    }
    
    // Make the kalman track fit and initialise with the seed track
    KalmanTrack* kalmanTrack = new KalmanTrack(track);
    double deltaChi2 = 0.;
    
    // Add each of the new clusters, and if the delta chi2 exceeds the cut then throw away this track extension
    KDCluster* kdStart = (*trackSegments[nTrackExtension])[0]->getEnd();
    double newHitDeltaChi2 = kalmanTrack->addCluster(kdStart);
    cout<<"- delta chi2 = "<<newHitDeltaChi2<<endl;
    if(newHitDeltaChi2 > 100.){
      delete kalmanTrack;
      delete trackSegments[nTrackExtension];
      continue;
    }
    
    deltaChi2+=newHitDeltaChi2;
    npoints++;
    
    for(unsigned int trackCell=0;trackCell<trackSegments[nTrackExtension]->size()-2;trackCell++){
      KDCluster* kdEnd = (*trackSegments[nTrackExtension])[trackCell]->getStart();
      newHitDeltaChi2 = kalmanTrack->addCluster(kdEnd);
      cout<<"- delta chi2 = "<<newHitDeltaChi2<<endl;

      if(newHitDeltaChi2 > 100.){
//        delete kalmanTrack;
        break;
      }

      deltaChi2+=newHitDeltaChi2;
      npoints++;

    }
    
    std::cout<<"-- final chi2/ndof increase: "<<deltaChi2/npoints<<" for "<<npoints<<" points"<<std::endl;
    if(bestTrack == NULL){
      bestTrack = kalmanTrack;
      bestChi2=deltaChi2;
      bestNpoints=npoints;
    }else if(deltaChi2/npoints < bestChi2/bestNpoints){
      delete bestTrack;
      bestTrack = kalmanTrack;
      bestChi2=deltaChi2;
      bestNpoints=npoints;
    }else{
      delete kalmanTrack;
    }
  }
  
  if(bestTrack != NULL){
    
    std::cout<<"== taking track with increase of chi2: "<<bestChi2<<", with "<<bestNpoints<<" points"<<std::endl;
    track->setKalmanTrack(bestTrack);
    
    for(int newHits=0; newHits < bestNpoints; newHits++){
      track->add( bestTrack->m_kalmanClusters[bestTrack->m_kalmanClusters.size()-1-newHits] );
    }
  }
  
  
  /*
//  std::cout<<"== Attempting to extend track"<<std::endl;
  // Get the inital track chi2/ndof
  int npoints = track->nPoints();
  double chi2ndof = track->chi2ndof();

  // Of all of the track segments, get the one with lowest chi2/ndof. Now look at each point on the cellular
  // track and add it to the track. If the delta chi2/ndof is small enough, keep the hit
  vector<KDTrack*> fittedTrackSegments;
  getFittedTracks(fittedTrackSegments, trackSegments, usedCells);
  KDTrack* bestTrackSegment = fittedTrackSegments[0];
  
  double bestChi2=1e10;
  KDTrack* bestExtension;
//  std::cout<<"- Have "<<fittedTrackSegments.size()<<" fitted segments. Original chi2/ndof is "<<chi2ndof<<std::endl;
  for(int nTrackExtension=0;nTrackExtension<fittedTrackSegments.size();nTrackExtension++){
    
    double newChi2ndof = fitWithExtension(*track, fittedTrackSegments[nTrackExtension]->m_clusters);
    
//    std::cout<<"Extension "<<nTrackExtension<<" gives a chi2/ndof of "<<newChi2ndof<<std::endl;
//    if( (newChi2ndof-bestChi2) < 0.5*bestChi2 ){
    if( newChi2ndof < bestChi2 ){
      bestChi2 = newChi2ndof;
      bestExtension = fittedTrackSegments[nTrackExtension];
    }

  }

//  if( bestChi2 != chi2ndof && ((bestChi2-chi2ndof) < 0.5*chi2ndof) ){
  if( bestChi2 != 1e10 && ((bestChi2-chi2ndof) < 10.*m_chi2cut) ){
    for(int newpoint=(bestExtension->m_clusters.size()-3);newpoint>=0;newpoint--){
      track->insert(bestExtension->m_clusters[newpoint]);
      if(bestChi2 < 10.) used[bestExtension->m_clusters[newpoint]] = true;
    }
  }
  
//  double newChi2ndof = fitWithExtension(*track, bestTrackSegment->m_clusters);

//  if( newChi2ndof < m_chi2cut ){
//    for(int newpoint=(bestTrackSegment->m_clusters.size()-3);newpoint>=0;newpoint--){
//      track->insert(bestTrackSegment->m_clusters[newpoint]);
//      if(newChi2ndof < 10.) used[bestTrackSegment->m_clusters[newpoint]] = true;
//    }
//  }
 //
}
*/

void ConformalTracking::extendTrack(KDTrack* track,std::vector<cellularTrack*> trackSegments, std::map<KDCluster*,bool>& used, std::map<Cell*,bool>& usedCells){
  
//  cout<<"== extending track, have "<<trackSegments.size()<<" candidates"<<endl;
  // For each track segment, perform a kalman filter on the addition points and chose the track extension with the
  // best delta chi2.
  double bestChi2=1.e9; double bestNpoints; vector<KDCluster*> finalHits;
  double bestChi2Conformal, bestChi2ZS;

//  for(int i=0;i<track->clusters().size();i++) std::cout<<"- cluster "<<i<<" has u,v "<<track->clusters()[i]->getU()<<","<<track->clusters()[i]->getV()<<std::endl;
  
  
  for(int nTrackExtension=0;nTrackExtension<trackSegments.size();nTrackExtension++){
    
//    std::cout<<"Track extension "<<nTrackExtension<<std::endl;
    vector<KDCluster*> newHits;
    
    // Add each of the new clusters
    KDCluster* kdStart = (*trackSegments[nTrackExtension])[0]->getEnd();
    newHits.push_back(kdStart);
//    std::cout<<"- cluster "<<" has u,v "<<kdStart->getU()<<","<<kdStart->getV()<<std::endl;
    
    for(unsigned int trackCell=0;trackCell<trackSegments[nTrackExtension]->size()-2;trackCell++){
      KDCluster* kdEnd = (*trackSegments[nTrackExtension])[trackCell]->getStart();
      newHits.push_back(kdEnd);
//      std::cout<<"- cluster "<<" has u,v "<<kdEnd->getU()<<","<<kdEnd->getV()<<std::endl;

    }
    
    double newChi2Conformal, newChi2ZS;
    double newChi2 = fitWithExtension(*track,newHits,newChi2Conformal,newChi2ZS);
    
    if(newChi2<bestChi2){
      finalHits.clear();
      for(int i=0;i<newHits.size();i++) finalHits.push_back(newHits[i]);
      bestChi2 = newChi2;
      bestChi2Conformal = newChi2Conformal;
      bestChi2ZS = newChi2ZS;
    }
    
  }
  if(finalHits.size() == 0) return;
  
//  std::cout<<"-- final chi2ZS/ndof increase: "<<(bestChi2ZS-track->chi2ndofZS())/finalHits.size()<<" for "<<finalHits.size()<<" points"<<std::endl;
//  std::cout<<"-- final chi2/ndof increase: "<<(bestChi2Conformal-track->chi2ndof())/finalHits.size()<<" for "<<finalHits.size()<<" points"<<std::endl;
  
  if((bestChi2ZS-track->chi2ndofZS()) < 100. && (bestChi2Conformal-track->chi2ndof()) < 100.){
    for(int i=finalHits.size()-1;i>=0;i--){
//      used[finalHits[i]] = true;
      track->m_clusters.insert(track->m_clusters.begin(),finalHits[i]);
      track->m_nPoints++;
    }
//    std::cout<<"- increasing chi2ZS/ndof from "<<track->chi2ndofZS()<<" to "<<bestChi2ZS<<std::endl;
//    std::cout<<"- increasing chi2/ndof from "<<track->chi2ndof()<<" to "<<bestChi2Conformal<<std::endl;
    track->linearRegression();
    track->linearRegressionConformal();
//    std::cout<<"- new chi2/ndof: "<<track->chi2ndof()<<", chi2ZS/ndof: "<<track->chi2ndofZS()<<std::endl;
  }
 
}

double ConformalTracking::fitWithExtension(KDTrack track, std::vector<KDCluster*> hits, double& newChi2, double& newChi2ZS){
  
  // Add the point to the track
  for(int i=0;i<hits.size();i++) track.add(hits[i]);

  int npoints = track.nPoints();
  /*
  // Fit the track and get the new chi2
      ROOT::Math::Functor FCNFunction(track,2);
      newFitter.SetFunction(FCNFunction);
//  globalTrack = &track;
  newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
  newFitter.SetVariable(1,"intercept", 0., 0.1);
  
  // Fit the track, first in uv space, then in sz space
  track.setConformalFit(true);
  newFitter.Minimize();
		
  // Now set the track parameters from the conformal fit, and fit in sz
  track.setGradient(newFitter.X()[0]);
  track.setIntercept(newFitter.X()[1]); */
//  track.setGradientError(newFitter.Errors()[0]);
//  track.setInterceptError(newFitter.Errors()[1]);
//  track.setConformalFit(false);
//  newFitter.Minimize();
//  track.setGradientZS(newFitter.X()[0]);
//  track.setInterceptZS(newFitter.X()[1]);
  
  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal();
//  track.fit();
  double chi2ndof = sqrt(track.chi2ndofZS()*track.chi2ndofZS() + track.chi2ndof()*track.chi2ndof());
  newChi2 = track.chi2ndof();
  newChi2ZS = track.chi2ndofZS();
  return chi2ndof;

}

double ConformalTracking::fitWithPoint(KalmanTrack kalTrack, KDCluster* hit){

//  KalmanTrack* kalTrack = new KalmanTrack(&track);
  
  double chi2ndof = kalTrack.addCluster(hit);
  
//  delete kalTrack;
  /*
  // Add the point to the track
  track.add(hit);
  
  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  if(track.chi2ndof() > 1000.) return track.chi2ndof();
  track.linearRegressionConformal();
  double chi2ndof = track.chi2ndof()+track.chi2ndofZS();
*/
  return chi2ndof;

}

double ConformalTracking::fitWithPoint(KDTrack kdTrack, KDCluster* hit){
  
  double chi2 = kdTrack.chi2ndof();
  double chi2zs = kdTrack.chi2ndofZS();
  kdTrack.add(hit);
  kdTrack.linearRegression();
  kdTrack.linearRegressionConformal(); //FCC study
  double newchi2 = kdTrack.chi2ndof();
  double newchi2zs = kdTrack.chi2ndofZS();

  return (newchi2-chi2 + newchi2zs-chi2zs);

}


// Debug function - checks if a track will be associated to an MC particle or not
double ConformalTracking::checkReal(KDTrack* track, std::map<KDCluster*,MCParticle*> kdParticles, std::map<MCParticle*,bool>& reconstructed, std::map<MCParticle*, std::vector<KDCluster*> > MCparticleHits){
 
  // Store all mcparticles associated to this track
  std::vector<MCParticle*> particles;
  std::map<MCParticle*,double> particleHits;
  double nHits=0.;
  
  // Get the clusters from this track
  std::vector<KDCluster*> clusters = track->m_clusters;
  
  // Loop over all hits and see which particle they are associated to
  for(int itCluster=0;itCluster<clusters.size();itCluster++){
    
    // Get the hit
    KDCluster* cluster = clusters[itCluster];
    nHits++;
    
    // If we already have hits on this particle, then just increment the counter
    if(particleHits.count(kdParticles[cluster]))
      particleHits[kdParticles[cluster]]++;
    else{
      // Otherwise register the new particle and set the number of hits to 1
      particles.push_back(kdParticles[cluster]);
      particleHits[kdParticles[cluster]] = 1;
    }
    
    
  }
  
  // Now look how many hits are on each particle and calculate the purity
  double bestHits = 0.; MCParticle* bestParticle=NULL;
  for(int iPart=0;iPart<particles.size();iPart++){
    if(particleHits[particles[iPart]] > bestHits){
      bestHits = particleHits[particles[iPart]];
      bestParticle = particles[iPart];
    }
  }
  
  // Calculate the purity
  double purity = bestHits/nHits;
  std::cout<<"Number of hits on track: "<<nHits<<". Good hits: "<<bestHits<<". Purity: "<<purity<<". Pt: "<<sqrt( bestParticle->getMomentum()[0]*bestParticle->getMomentum()[0] + bestParticle->getMomentum()[1]*bestParticle->getMomentum()[1] )<<". Track chi2/ndof: "<<track->chi2ndof()<<". Chi2/ndof in SZ fit: "<<track->chi2ndofZS()<<std::endl;
  
  // Check if any hits are missing
  std::vector<KDCluster*> mcHits = MCparticleHits[bestParticle];
  int uniqueHits = getUniqueHits(mcHits);
  
  std::cout<<"Track should contain "<<uniqueHits<<" hits, "<< ((uniqueHits > bestHits ) ? "is missing hits" : "all hits found.") << std::endl;

  
  std::vector<KDCluster*> trackHits = track->m_clusters;
  for(int th=0;th<trackHits.size();th++) std::cout<<"Hit "<<th<<" u = "<<trackHits[th]->getU()<<" v = "<<trackHits[th]->getV()<<" x = "<<trackHits[th]->getX()<<" y = "<<trackHits[th]->getY()<<" z = "<<trackHits[th]->getZ()<<std::endl;

//  for(int itCluster=0;itCluster<clusters.size();itCluster++) std::cout<<"Hit "<<itCluster<<" has position: "<<clusters[itCluster]->getU()<<","<<clusters[itCluster]->getV()<<std::endl;
  std::cout<<"== Terms in conformal fit: "<<track->m_intercept<<", "<<track->m_gradient<<", "<<track->m_quadratic<<std::endl;
  std::cout<<"== Terms in zs fit: "<<track->m_interceptZS<<", "<<track->m_gradientZS<<std::endl;

  track->linearRegressionConformal(true);
  track->calculateChi2SZ(NULL,true);
  if(purity >= m_purity){
    reconstructed[bestParticle] = true;
  }
  
  // Return the purity
  return purity;
  
}

// Debug function - gets number of unique hits
int ConformalTracking::getUniqueHits(std::vector<KDCluster*> hits){
  
  int nUniqueHits = 0;
  
  for(int iHit=0;iHit<hits.size();iHit++){
   
    // Get the cluster
    KDCluster* hit = hits[iHit];
    bool sameID = false;

    // Check if any following clusters have the same detector id
    for(int iHitLater=iHit+1;iHitLater<hits.size();iHitLater++){
      
      // Get the new cluster
      KDCluster* nhit = hits[iHitLater];

      if(hit->sameLayer(nhit)) sameID = true;
    }
    
    // Didn't find another hit with the same ID
    if(!sameID) nUniqueHits++;
    
  }
  
  return nUniqueHits;
}

//---- Make the track in the same way as the pattern recognition would, but
//---- without the inclusion of unassociated hits. See which criteria fail

void ConformalTracking::checkReconstructionFailure(MCParticle* particle, std::map<MCParticle*, std::vector<KDCluster*> > particleHits, std::map<KDCluster*,bool> used, KDTree* nearestNeighbours){
  
  
  // Get the hits for this MC particle
  std::vector<KDCluster*> trackHits = particleHits[particle];

  // Sort the hits from larger to smaller radius
  std::sort(trackHits.begin(),trackHits.end(),sort_by_radiusKD);
  std::cout<<"Track contains "<<getUniqueHits(trackHits)<<" unique hits. Track size is "<<trackHits.size()<<std::endl;

    for(int th=0;th<trackHits.size();th++) std::cout<<"Hit "<<th<<" u = "<<trackHits[th]->getU()<<" v = "<<trackHits[th]->getV()<<" x = "<<trackHits[th]->getX()<<" y = "<<trackHits[th]->getY()<<" z = "<<trackHits[th]->getZ()<<". "<<(trackHits[th]->used() ? "Used!" : "")<<std::endl;

  // Take the seed hit and build cell to 2nd
  KDCluster* seedHit = trackHits[0];
  VecCluster results;
  std::cout<<"- getting nearest neighbours for seed"<<std::endl;
//  nearestNeighbours->allNeighboursInRadius(seedHit, m_maxDistance, results);
  double theta = seedHit->getTheta();
  nearestNeighbours->allNeighboursInTheta(theta, m_thetaRange, results);
  std::cout<<"- got nearest neighbours for seed"<<std::endl;
  
  // If we don't produce the seed cell we have failed
  if( std::find(results.begin(),results.end(),trackHits[1]) == results.end() ){
    std::cout<<"- seed cell not produced, neighbour not inside search window"<<std::endl;
    double deltaU = fabs(trackHits[1]->getU()-trackHits[0]->getU());
    double deltaV = fabs(trackHits[1]->getV()-trackHits[0]->getV());
    std::cout<<"- distance between hits is "<< sqrt(deltaU*deltaU + deltaV*deltaV)<<" and search distance is "<<m_maxDistance<<std::endl;
    std::cout<<"- theta of hits is "<<trackHits[0]->getTheta()<<" and "<<trackHits[1]->getTheta()<<", delta theta = "<<trackHits[0]->getTheta()-trackHits[1]->getTheta()<<std::endl;
//    return;
  }

  if(trackHits[1]->getR() >= trackHits[0]->getR()){
    std::cout<<"- seed cell not produced, neighbour at higher radius"<<std::endl;
//    return;
  }
  
  double length = sqrt( (trackHits[1]->getU()-trackHits[0]->getU())*(trackHits[1]->getU()-trackHits[0]->getU()) + (trackHits[1]->getV()-trackHits[0]->getV())*(trackHits[1]->getV()-trackHits[0]->getV()) );
  if(length > m_maxDistance){
    std::cout<<"- seed cell not produced, neighbour too far away"<<std::endl;
//    return;
  }
  
  
  // Make the seed cell and extrapolate until all hits found
  Cell* seedCell = new Cell(trackHits[0],trackHits[1]);
  std::vector<Cell*> cells; cells.push_back(seedCell);
  std::cout<<"have seed cell"<<std::endl;
  // CHECK CELL LENGTH!! TO DO

  int hitNumber = 1; // Current hit after the seed cell has been formed
  bool trackBuilding = true;
  while(hitNumber<(trackHits.size()-1)){
    
    // Extend the current cell
    double searchDistance = m_maxDistance; //hit->getR();
//    if(searchDistance > trackHits[hitNumber]->getR()) searchDistance = trackHits[hitNumber]->getR();
    
    // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
    KDCluster* fakeHit = extrapolateCell(cells.back(),searchDistance/2.); // TODO: make this search a function of radius
    VecCluster results2;
    nearestNeighbours->allNeighboursInRadius(fakeHit, 1.25*searchDistance/2., results2);
    delete fakeHit;
    
    // Check if we found the hit we wanted
    if( std::find(results2.begin(),results2.end(),trackHits[hitNumber+1]) == results2.end() ){
      std::cout<<"- could not find hit number "<<hitNumber+1<<" inside the search window"<<std::endl;
      double deltaU = fabs(trackHits[hitNumber+1]->getU()-trackHits[hitNumber]->getU());
      double deltaV = fabs(trackHits[hitNumber+1]->getV()-trackHits[hitNumber]->getV());
      std::cout<<"- distance between hits is "<< sqrt(deltaU*deltaU + deltaV*deltaV)<<" and search distance is "<<m_maxDistance<<std::endl;
//      for(int c=0;c<cells.size();c++) delete cells[c];
//      return;
    }
    
    // Check radial conditions
    if(trackHits[hitNumber+1]->getR() >= trackHits[hitNumber]->getR()){
      std::cout<<"- cell "<<hitNumber<<" not produced, neighbour at higher radius"<<std::endl;
//      return;
    }
    
    // Make the cell for extrapolation in the next round
    Cell* cell = new Cell(trackHits[hitNumber],trackHits[hitNumber+1]);
    std::cout<<"have new cell"<<std::endl;

    // Check if the cell would be killed by cuts
    if( cells.back()->getAngle(cell) > m_maxCellAngle || cells.back()->getAngleRZ(cell) > m_maxCellAngleRZ ){
      if(cells.back()->getAngle(cell) > m_maxCellAngle) std::cout<<"- cell "<<hitNumber<<" killed by angular cut. Angle is "<<cells.back()->getAngle(cell)<<std::endl;
      if(cells.back()->getAngleRZ(cell) > m_maxCellAngleRZ) std::cout<<"- cell "<<hitNumber<<" killed by RZ angular cut. Angle is "<<cells.back()->getAngleRZ(cell)<<std::endl;
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
  std::map<Cell*,bool> usedCells;
  createTracksNew(cellularTracks, cells.back(), usedCells);

  if(cellularTracks.size() != 1){
    std::cout<<"- cellular track not produced from cells. Returned "<<cellularTracks.size()<<" candidates"<<std::endl;
    for(int c=0;c<cells.size();c++) delete cells[c];
    for(int ct=0;ct<cellularTracks.size();ct++) delete cellularTracks[ct];
    return;
  }

  // Check that all of the hits are on the track
  if(cellularTracks[0]->size() != (trackHits.size()-1)){
    std::cout<<"- failed to put all cells on track. Cellular track size: "<<cellularTracks[0]->size()<<std::endl;
    std::cout<<"- number of cells held: "<<cells.size()<<std::endl;
    return;
  }

  // Now test the built-in functions to make KDTracks from CellularTracks
  std::vector<KDTrack*> finalTracks;
  getFittedTracks(finalTracks,cellularTracks,usedCells);

  if(finalTracks.size() != 1){
    std::cout<<"- kd track not produced from cellular track. Returned "<<finalTracks.size()<<" candidates"<<std::endl;
    for(int c=0;c<cells.size();c++) delete cells[c];
    for(int t=0;t<finalTracks.size();t++) delete finalTracks[t];
    return;
  }

  // Test the chi2 criteria
  KDTrack* mcTrack = finalTracks[0];
    mcTrack->linearRegressionConformal(true);

  std::cout<<"== Track chi2/ndof is "<<mcTrack->chi2ndof()<<", ZS chi2/ndof is "<<mcTrack->chi2ndofZS()<<std::endl;
  std::cout<<"== Terms in conformal fit: "<<mcTrack->m_intercept<<", "<<mcTrack->m_gradient<<", "<<mcTrack->m_quadratic<<std::endl;
  std::cout<<"== Terms in zs fit: "<<mcTrack->m_interceptZS<<", "<<mcTrack->m_gradientZS<<std::endl;
  mcTrack->calculateChi2SZ(NULL,true);

  if(mcTrack->chi2ndof() > m_chi2cut || mcTrack->chi2ndofZS() > m_chi2cut) {
    if(mcTrack->chi2ndof() > m_chi2cut) std::cout<<"- track killed by chi2/ndof cut. Track chi2/ndof is "<<mcTrack->chi2ndof()<<std::endl;
    if(mcTrack->chi2ndofZS() > m_chi2cut) std::cout<<"- track killed by ZS chi2/ndof cut. Track chi2/ndof in ZS is "<<mcTrack->chi2ndofZS()<<std::endl;
    for(int c=0;c<cells.size();c++) delete cells[c];
    for(int t=0;t<finalTracks.size();t++) delete finalTracks[t];
    return;
  }

  std::cout<<"== Should have produced this track"<<std::endl;
  for(int c=0;c<cells.size();c++) delete cells[c];
  for(int t=0;t<finalTracks.size();t++) delete finalTracks[t];
  return;

}







