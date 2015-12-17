/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
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
//TGraphErrors* fitter;

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
  registerProcessorParameter( "DebugPlots", 				"Plots for debugging the tracking", 													m_debugPlots, 				bool(false) 	);
  registerProcessorParameter( "ThetaRange", 				"Angular range for initial cell seeding", 										m_thetaRange, 				double(0.1)		);
  registerProcessorParameter( "MaxCellAngle", 			"Cut on angle between two cells for cell to be valid", 				m_maxCellAngle, 			double(0.035)	);
  registerProcessorParameter( "MaxDistance",				"Maximum length of a cell (max. distance between two hits)", 	m_maxDistance, 				double(0.015) 	);
  registerProcessorParameter( "MaxChi2",						"Maximum chi2/ndof for linear conformal tracks", 							m_chi2cut, 						double(300.) 	);
  registerProcessorParameter( "MinClustersOnTrack", "Minimum number of clusters to create a track", 							m_minClustersOnTrack, int(5)				);
  
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

// Sort cells from higher to lower weight
bool sort_by_cellWeight(Cell* cell1, Cell* cell2){
  int weight1 = cell1->getWeight();
  int weight2 = cell2->getWeight();
  return (weight1 > weight2);
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
  
  // Initialise histograms (if debug plotting on)
  if(m_debugPlots){
    
    // Automatically add histograms to output root file
    AIDAProcessor::histogramFactory(this);
    
    // Histogram initailisation
    
    // Histograms for tuning parameters (cell angle cut, cell length cut)
    m_cellAngle = new TH1F("cellAngle","cellAngle",1250,0,0.05);
    m_cellAngleRadius = new TH2F("cellAngleRadius","cellAngleRadius",400,0,0.04,1000,0,0.04);
    m_cellLengthRadius = new TH2F("cellLengthRadius","cellLengthRadius",300,0,0.03,1000,0,0.04);
    m_cellAngleLength = new TH2F("cellAngleLength","cellAngleLength",400,0,0.04,300,0,0.03);
    m_conformalChi2 = new TH1F("conformalChi2","conformalChi2",100,0,100);
    m_conformalChi2real = new TH1F("conformalChi2real","conformalChi2real",1000,0,1000);
    m_conformalChi2fake = new TH1F("conformalChi2fake","conformalChi2fake",1000,0,1000);
    m_conformalChi2Purity = new TH2F("conformalChi2Purity","conformalChi2Purity",150,0,1.5,1000,0,1000);
    
    m_cellAngleMC = new TH1F("cellAngleMC","cellAngleMC",1250,0,0.05);
    m_cellAngleRadiusMC = new TH2F("cellAngleRadiusMC","cellAngleRadiusMC",400,0,0.04,1000,0,0.04);
    m_cellLengthRadiusMC = new TH2F("cellLengthRadiusMC","cellLengthRadiusMC",300,0,0.03,1000,0,0.04);
    m_cellAngleLengthMC = new TH2F("cellAngleLengthMC","cellAngleLengthMC",400,0,0.04,300,0,0.03);
    
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
  
	// Set up the conformal fitter
  newFitter.SetMaxFunctionCalls(1000000);
  newFitter.SetMaxIterations(100000);
  newFitter.SetPrecision(1e-10);
  newFitter.SetPrintLevel(0);
  gErrorIgnoreLevel=kWarning;
  
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
  UTIL::BitField64 m_encoder( lcio::ILDCellID0::encoder_string ) ;

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
    getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;
    
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
        // Make the conformal hit
        KDCluster* kdhit = new KDCluster(hit);
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
//      if(particlePt < 0.1) continue;
      // Cut on stable particles
      if(mcParticle->getGeneratorStatus() != 1) continue;
      // Sort the hits from larger to smaller radius
      std::sort(trackHits.begin(),trackHits.end(),sort_by_radiusKD);
      
//      if(particlePt>29.5 && particlePt<29.8){
//        debugHits = trackHits;
//        std::cout<<"==> Picked up "<<debugHits.size()<<" DEBUG hits"<<std::endl;
//      }
      // Loop over all hits for debugging
      for(int itHit=0;itHit<trackHits.size();itHit++){
        
        // Get the conformal clusters
        KDCluster* cluster = trackHits[itHit];
        
//        if(particlePt>29.5 && particlePt<29.8){
//          debugHits = trackHits;
//          std::cout<<"==> Picked up "<<debugHits.size()<<" DEBUG hits"<<std::endl;
//          debugHitCollection->addElement(tempHolder[cluster]);
//        }
        
        
//        const int celId = tempHolder[cluster]->getCellID0() ;
//        m_encoder.setValue(celId) ;
//        int subdet = m_encoder[lcio::ILDCellID0::subdet];
//        int side = m_encoder[lcio::ILDCellID0::side];
//        int layer = m_encoder[lcio::ILDCellID0::layer];
//        
//        double x = tempHolder[cluster]->getPosition()[0];
//        double y = tempHolder[cluster]->getPosition()[1];
//        double z = tempHolder[cluster]->getPosition()[2];
//        double r = sqrt(x*x + y*y);
        //        std::cout<<"Subdetector "<<subdet<<", side "<<side<<", layer "<<layer<<std::endl;
        //        std::cout<<"Hit position "<<x<<", "<<y<<", "<<z<<". radius "<<r<<std::endl;
        
      }
      
      
      // Now loop over the hits and make cells - filling histograms along the way
      int nHits = trackHits.size();
      //      if(nHits > 10){std::cout<<"Skipping particle"<<std::endl; continue;}
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
        double cell0Length = sqrt( pow(cluster0->getU() - cluster1->getU(),2) + pow(cluster0->getV() - cluster1->getV(),2) );
        double cell1Length = sqrt( pow(cluster1->getU() - cluster2->getU(),2) + pow(cluster1->getV() - cluster2->getV(),2) );
        
        m_cellAngleMC->Fill(angleBetweenCells);
        m_cellAngleRadiusMC->Fill(cluster2->getR(),angleBetweenCells);
        m_cellLengthRadiusMC->Fill(cluster0->getR(),cell0Length);
        m_cellAngleLengthMC->Fill(cell1Length,angleBetweenCells);
        
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
            //            if(cell0Length > 0.02) cout<<"Cell 0 is very long!"<<endl;
            drawline(cluster0,cluster1,itHit+1);
          }
          //          if(cell1Length > 0.02) cout<<"Cell "<<itHit+1<<" is very long!"<<endl;
          // Draw line style differently if the cell angle was too large
          //          if(angleBetweenCells > (m_maxCellAngle*exp(-0.0015/cluster2->getR()))){
          if(angleBetweenCells > (m_maxCellAngle)){
            drawline(cluster1,cluster2,itHit+2,3);
//            std::cout<<"Cell angle of "<<angleBetweenCells<<" is too large"<<std::endl;
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
  
//  std::vector<KDCluster*> kdClusters; 								// Conformal hits
  std::map<KDCluster*,TrackerHitPlane*> kdClusterMap; 	// Their link to "real" hits
  std::map<KDCluster*,bool> used;												// Map of whether a hit has been included in a track or not
  std::map<KDCluster*,bool> used2;
  std::vector<KDTrack> conformalTracks;									// KD tracks - each is a list of kd hits in the found tracks
  std::vector<double> conformalTracksChi2ndof;					// A list of chi2/ndof for all conformal tracks

  // Loop over all input collections. Tracking will be attempted on the collection, then hits from the next collection
  // will be added to the unused hits already there.
  std::map< int,std::vector<KDCluster*> > collectionClusters;
  for(unsigned int collection=0; collection<trackerHitCollections.size();collection++){
    
    // Loop over tracker hits and make conformal hit collection
    std::vector<KDCluster*> tempClusters;
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    std::cout<<"Collection "<<collection<<" has "<<nHits<<" hits"<<std::endl;
    for(int itHit=0;itHit<nHits;itHit++){
      
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
      
      // Make a new kd cluster (if debugging then these already exist - don't remake them)
      KDCluster* kdhit;
      if(m_debugPlots) kdhit = conformalHits[hit];
      if(!m_debugPlots) kdhit = new KDCluster(hit);
      
      // Get subdetector information and add it to the kdhit
      const int celId = hit->getCellID0() ;
      m_encoder.setValue(celId) ;
      int subdet = m_encoder[lcio::ILDCellID0::subdet];
      int side = m_encoder[lcio::ILDCellID0::side];
      int layer = m_encoder[lcio::ILDCellID0::layer];
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
  
  for(unsigned int collection=0; collection<trackerHitCollections.size();collection++){
    
    /*
    // Loop over tracker hits and make conformal hit collection
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    for(int itHit=0;itHit<nHits;itHit++){
      
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
      
      // Make a new kd cluster (if debugging then these already exist - don't remake them)
      KDCluster* kdhit;
      if(m_debugPlots) kdhit = conformalHits[hit];
      if(!m_debugPlots) kdhit = new KDCluster(hit);
      
      // Get subdetector information and add it to the kdhit
      const int celId = hit->getCellID0() ;
      m_encoder.setValue(celId) ;
      int subdet = m_encoder[lcio::ILDCellID0::subdet];
      int side = m_encoder[lcio::ILDCellID0::side];
      int layer = m_encoder[lcio::ILDCellID0::layer];
      kdhit->setDetectorInfo(subdet,side,layer);
      
      // Store the link between the two
      kdClusterMap[kdhit] = hit;
      kdClusters.push_back(kdhit);
      
      // Debug histogramming
      if(m_debugPlots && m_eventNumber == 0){
        m_conformalEvents->Fill( kdhit->getU(),kdhit->getV() );
        m_nonconformalEvents->Fill( hit->getPosition()[0], hit->getPosition()[1] );
        m_conformalEventsRTheta->Fill( kdhit->getR(), kdhit->getTheta() );
      }
      
    }
    */
    std::vector<KDCluster*> kdClusters; 								// Conformal hits

    for(unsigned int col=0;col<=collection;col++){
      
      std::vector<KDCluster*> clusters = collectionClusters[col];
      int nhits = clusters.size();
      for(int hit=0;hit<nhits;hit++){
        if(used.count(clusters[hit])) continue;
        kdClusters.push_back(clusters[hit]);
      }
      
    }
    
    // Sort the KDClusters from larger to smaller radius
    if(kdClusters.size() == 0) continue;
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
    
    // Loop over all current tracks
    for(int currentTrack=0;currentTrack<nCurrentTracks;currentTrack++){

      // This step (although first) only runs when tracks have already been produced. An attempt
      // is made to extend them with hits from the new collection, using the final cell of the track
      // as the seed cell.
      
      // Containers to hold new cells made, and to check if a hit already has a cell connected to it
      std::vector<Cell*> cells;
      
      // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
      Cell* seedCell = new Cell(conformalTracks[currentTrack].clusters()[1],conformalTracks[currentTrack].clusters()[0]);
      cells.push_back(seedCell);
      
      // All seed cells have been created, now try create all "downstream" cells until no more can be added
      extendSeedCells(cells, used, nearestNeighbours, false, debugHits);
      
      // We create all acceptable tracks by looping over all cells with high enough weight to create
      // a track and trace their route back to the seed hit. We then have to choose the best candidate
      // at the end (by minimum chi2 of a linear fit)
      std::map<Cell*,bool> usedCells;
      std::map<Cell*,bool> usedCells2;
      std::vector<cellularTrack> trackSegments;
      
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
        if(cells[itCell]->getWeight() < 2) break;
        
        // Produce all segments leading back to the track from this cell
        std::vector<cellularTrack> candidateSegments = createTracks(cells[itCell],usedCells2);
        
        // Store all of these segments for later
        if(candidateSegments.size() == 0) continue;
        trackSegments.insert(trackSegments.end(),candidateSegments.begin(), candidateSegments.end());

        // Mark the cells from these segments as having been used
        for(unsigned int itSegment=0;itSegment<candidateSegments.size();itSegment++){
          for(unsigned int itCell=0;itCell<candidateSegments[itSegment].size();itCell++){
            usedCells[candidateSegments[itSegment][itCell]]=true;
          }
        }
      }

      // Decide which segment to add on to the track, and mark the added hits as used
      if(trackSegments.size() == 0) continue;
      extendTrack(conformalTracks[currentTrack],trackSegments,used);
      
      // Clean up
      for(unsigned int itCell=0;itCell<cells.size();itCell++) delete cells[itCell];

    }
 
  	// ---------------------------------------------------------------------
		// Try to create new tracks using all of the kdHits currently held
	  // ---------------------------------------------------------------------
		unsigned int nKDHits = kdClusters.size();
    streamlog_out( DEBUG4 )<<"Seeding with hits"<<std::endl;
    streamlog_out( DEBUG4 )<<"Attempting to seed with hits: "<<nKDHits<<std::endl;
    
    // Loop over all current hits
    for(unsigned int nKDHit = 0; nKDHit<nKDHits; nKDHit++){
      
      // Get the kdHit and check if it has already been used (assigned to a track)
      KDCluster* kdhit = kdClusters[nKDHit];
      if(used.count(kdhit)) continue;
      if(kdhit->getR() < 0.005) break; // new cut - once we get to inner radius we will never make tracks. temp? TODO: make parameter?

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
      
      // Sort the neighbours from outer to inner radius
      if(results.size() == 0) continue;
      std::sort(results.begin(),results.end(),sort_by_radiusKD);
      
      // Objects to hold cells
      std::vector<Cell*> cells;
      
      // Make seed cells pointing inwards (conformal space)
      for(unsigned int neighbour=0;neighbour<results.size();neighbour++){
        
        // Get the neighbouring hit
        KDCluster* nhit = results[neighbour];
        
        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
        if(used.count(nhit)) continue;
        if(kdhit->sameLayer(nhit)) continue;
        if(nhit->getR() >= kdhit->getR()) continue;

        // If close to the central region increasing r != increasing z
        double zdifference = abs(kdhit->getZ()-nhit->getZ());
        if(zdifference > 10.){
          if(kdhit->getZ() > 0. && nhit->getZ() < kdhit->getZ()) continue;
          if(kdhit->getZ() < 0. && nhit->getZ() > kdhit->getZ())continue;
        }

        // Create the new seed cell
        Cell* cell = new Cell(kdhit,nhit);
        cells.push_back(cell);
        
        // Debug plotting
        if(m_debugPlots && m_eventNumber == 0){
          m_canvConformalEventDisplayAllCells->cd();
          drawline(kdhit,nhit,1);
        }

      }
      
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
      std::vector<KDTrack> cellTracks;
      std::vector<double> cellTracksChi2ndof;
      
      // Sort Cells from highest to lowest weight
      std::sort(cells.begin(),cells.end(),sort_by_cellWeight);

      // Create tracks by following a path along cells
      int nCells = cells.size();
      for(int itCell=0;itCell<nCells;itCell++){
        
        // Check if this cell has already been used
        if(usedCells.count(cells[itCell])) continue;

        // Check if this cell could produce a track (is on a long enough chain)
        if(cells[itCell]->getWeight() < (m_minClustersOnTrack-2)) break;
        
        // Produce all tracks leading back to the seed hit from this cell
        std::vector<cellularTrack> candidateTracks = createTracks(cells[itCell],usedCells2);

        // Debug plotting
        if(m_debugPlots && m_eventNumber == 0){
          m_canvConformalEventDisplayAcceptedCells->cd();
          for(int iTr=0;iTr<candidateTracks.size();iTr++){
            cellularTrack track = candidateTracks[iTr];
            for(unsigned int trackCell=0;trackCell<track.size();trackCell++){
              drawline(track[trackCell]->getStart(),track[trackCell]->getEnd(),track.size()-trackCell);
            }
          }
        }

        // Look at the candidate tracks and fit them (+pick the best chi2/ndof)
        if(candidateTracks.size() == 0) continue;
        std::vector<double> chi2ndof;
        std::vector<KDTrack> bestTracks = getFittedTracks(candidateTracks,chi2ndof);

        // Store track(s) for later
        cellTracks.insert(cellTracks.end(),bestTracks.begin(),bestTracks.end());
        cellTracksChi2ndof.insert(cellTracksChi2ndof.end(),chi2ndof.begin(), chi2ndof.end());
        
        // Mark the cells as having been used, they will not be re-used
        for(unsigned int acceptedCandidates=0;acceptedCandidates<candidateTracks.size();acceptedCandidates++){
//          if(chi2ndof[acceptedCandidates] > 10.) continue;
        	for(unsigned int trackCell=0;trackCell<candidateTracks[acceptedCandidates].size();trackCell++){
          	usedCells[candidateTracks[acceptedCandidates][trackCell]]=true;
        	}
        }
        
      }
      
      // All tracks leading back to the seed hit have now been found. Decide which are the feasible candidates (may be more than 1)
      if(cellTracks.size() == 0) continue;
      std::vector<double> chi2ndof;
      std::vector<KDTrack> bestTracks = getLowestChi2(cellTracks,cellTracksChi2ndof,chi2ndof);
      
      // Could now think to do the full helix fit and apply a chi2 cut. TODO
  
      // Store the final CA tracks
      for(unsigned int itTrack=0;itTrack<bestTracks.size();itTrack++){
        
        // Cut on chi2
        if(chi2ndof[itTrack] > m_chi2cut) continue;
        
        // Check if the new track is a clone
        bool clone = false;
        for(int existingTrack=0;existingTrack<conformalTracks.size();existingTrack++){
          // If the track is a clone with lower chi2/ndof, replace the old track with the new one
          if( sameTrack(bestTracks[itTrack],conformalTracks[existingTrack]) ){
            clone = true;
            if( chi2ndof[itTrack] < conformalTracksChi2ndof[existingTrack] ){
              conformalTracks[existingTrack] = bestTracks[itTrack];
              conformalTracksChi2ndof[existingTrack] = chi2ndof[itTrack];
            }
          }
        }
        // If not a clone, save the new track
        if(!clone){
        	conformalTracks.push_back(bestTracks[itTrack]);
        	conformalTracksChi2ndof.push_back(chi2ndof[itTrack]);
        	std::cout<<"Pushing back best track with chi2/ndof "<<chi2ndof[itTrack]<<std::endl;
        }
        
        // If the track has a good chi2/ndof (so we are 'sure' that it is real) then mark the hits as used so that they are not used again
        if(chi2ndof[itTrack] < 10.){
          for(unsigned int itHit=0;itHit<bestTracks[itTrack].clusters().size();itHit++){
            used[bestTracks[itTrack].clusters()[itHit]] = true;
          }
        }
      }
      
      // Clean up
      for(unsigned int itCell=0;itCell<cells.size();itCell++) delete cells[itCell];

      // Fill debug plots
 /*     if(m_debugPlots){
        
        for(int itrack=0; itrack<cellTracks.size();itrack++){
          cellularTrack debugTrack = cellTracks[itrack];
          
          m_conformalChi2->Fill(chi2ndof);
          double purity = checkReal(debugTrack,kdParticles,reconstructed);
          if(purity >= 0.75){
            m_conformalChi2real->Fill(chi2ndof);
          }
          if(purity < 0.75){
            m_conformalChi2fake->Fill(chi2ndof);
          }
          m_conformalChi2Purity->Fill(purity,chi2ndof);
        }
      }*/
      
    }
//    }
    // Now finished looking at this hit collection. All tracks which can have extra hits added
    // have had them added, and no new tracks were possible using the sum of all collections
    // till now. Add the next collection and try again...
    delete nearestNeighbours;
  }
  
  // Now make "real" tracks from all of the conformal tracks
  std::cout<<"*** CA has made "<<conformalTracks.size()<<" tracks ***"<<std::endl;

  // Loop over all track candidates
  for(unsigned int caTrack=0;caTrack<conformalTracks.size();caTrack++){
    
    // Vector of all the hits on the track
    KDTrack conformalTrack = conformalTracks[caTrack];
    streamlog_out( DEBUG5 )<<"Made a track with "<<conformalTrack.clusters().size()<<" hits"<<std::endl;
    
    // Make the LCIO track hit vector
    EVENT::TrackerHitVec trackHits;
    for(unsigned int itHit=0;itHit<conformalTrack.clusters().size();itHit++){
      KDCluster* cluster = conformalTrack.clusters()[itHit];
      trackHits.push_back(kdClusterMap[cluster]);
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
//    }// delete track; delete marlinTrack; continue;}
    
    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ] = trackHits.size();
    
    // Push back to the output container
    trackCollection->addElement(track);
    
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
    
    int nReconstructed = 0;
    // Additionally draw all tracks that were not reconstructed
    m_canvConformalEventDisplayMCunreconstructed->cd();
    int nParticles = particleCollection->getNumberOfElements();
    for(int itP=0;itP<nParticles;itP++){
      // Get the particle
      MCParticle* mcParticle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itP) ) ;
      // Check if it was reconstructed
      if(reconstructed.count(mcParticle)){nReconstructed++; continue;}
      // Check if it was stable
      if(mcParticle->getGeneratorStatus() != 1) continue;
      // Get the conformal hits
      if(particleHits.count(mcParticle) == 0) continue;
      std::vector<KDCluster*> mcHits = particleHits[mcParticle];
      // Cut on the number of hits
      int uniqueHits = getUniqueHits(mcHits);
      if(uniqueHits < m_minClustersOnTrack) continue;
      std::sort(mcHits.begin(),mcHits.end(),sort_by_radiusKD);
      // Draw the cells connecting the hits
      for(int iHit=0;iHit<(mcHits.size()-1);iHit++){
        drawline(mcHits[iHit],mcHits[iHit+1],iHit+1);
      }
      // List the pt of unreconstructed particles
      double particlePt = sqrt( mcParticle->getMomentum()[0]*mcParticle->getMomentum()[0] + mcParticle->getMomentum()[1]*mcParticle->getMomentum()[1] );
//      std::cout<<"Unreconstructed particle pt: "<<particlePt<<std::endl;
    }
//    std::cout<<"Reconstructed "<<nReconstructed<<" particles in total"<<std::endl;
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
      if(searchDistance > hit->getR()) searchDistance = hit->getR();
      
      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      KDCluster* fakeHit = extrapolateCell(cells[itCell],searchDistance/2.); // TODO: make this search a function of radius
      VecCluster results;
      nearestNeighbours->allNeighboursInRadius(fakeHit, searchDistance/2., results);
      delete fakeHit;
      
      if(extendingTrack) std::cout<<"Found "<<results.size()<<" neighbours from cell extrapolation"<<std::endl;
      // Make new cells pointing inwards
      for(unsigned int neighbour=0;neighbour<results.size();neighbour++){
        
        if(extendingTrack) std::cout<<"looking at neighbour "<<neighbour<<std::endl;
        // Get the neighbouring hit
        KDCluster* nhit = results[neighbour];
        
        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away from IP
        if(used.count(nhit)){if(extendingTrack) std::cout<<"- used"<<std::endl; continue;}
        if(hit->sameLayer(nhit)){if(extendingTrack) std::cout<<"- same layer"<<std::endl; continue;}
        if(nhit->getR() >= hit->getR()){if(extendingTrack) std::cout<<"- higher radius"<<std::endl; continue;}
        double zdifference = abs(hit->getZ()-nhit->getZ());
        if(zdifference > 10.){
          if(hit->getZ() > 0. && nhit->getZ() < hit->getZ()){if(extendingTrack) std::cout<<"- z cut"<<std::endl; continue;}
          if(hit->getZ() < 0. && nhit->getZ() > hit->getZ()){if(extendingTrack) std::cout<<"- z cut"<<std::endl; continue;}
        }
        
        // Check if this cell already exists (rejoining branch)
        if(existingCells.count(hit) != 0){
          bool alreadyExists=false;
          int nExistingCells = existingCells[hit].size();
          for(int iterExisting=0;iterExisting<nExistingCells;iterExisting++){
            if( existingCells[hit][iterExisting]->getEnd() == nhit ){
              alreadyExists = true;
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
        if( cells[itCell]->getAngle(cell) > (m_maxCellAngle*exp(-0.001/nhit->getR())) ){
//        if( cells[itCell]->getAngle(cell) > m_maxCellAngle ){
          
          // Debug plotting
          if(m_debugPlots && m_eventNumber == 0){
            m_canvConformalEventDisplayAllCells->cd();
            drawline(hit,nhit,cells[itCell]->getWeight()+2,3);
          }

          delete cell;
          continue;
        }
        
        // Set the information about which cell this new cell is attached to and store it
        cell->setFrom(cells[itCell]);
        cells[itCell]->setTo(cell);
        cells.push_back(cell);
        existingCells[hit].push_back(cell);
        
        // Debug plotting
        if(m_debugPlots && m_eventNumber == 0){
          m_canvConformalEventDisplayAllCells->cd();
          drawline(hit,nhit,cells[itCell]->getWeight()+2);
        }
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

// TODO: check the logic in used cell notation, branching, weights of cells being followed. In particular if we increase weight and make continuously connected
// cells with non-continuous weight...
std::vector<cellularTrack> ConformalTracking::createTracks(Cell* seedCell, std::map<Cell*,bool>& usedCells){
  
  cellularTrack seedTrack;
  seedTrack.push_back(seedCell);
  
  std::vector<cellularTrack> cellularTracks;
  cellularTracks.push_back(seedTrack);
  
  while(toBeUpdated(cellularTracks)){
    int removeTrack=-1;
    int nTracks=cellularTracks.size();
    //    std::cout<<" -- check for update on "<<nTracks<<" tracks"<<std::endl;
    for(int iTrack=0;iTrack<nTracks;iTrack++){
      if(removeTrack > (-1)) continue;
      if( cellularTracks[iTrack].back()->getWeight() > 0 ) followPath(cellularTracks,iTrack,usedCells,removeTrack);
    }
    if(removeTrack > (-1)){
      //      std::cout<<"Erasing track "<<removeTrack<<std::endl;
      cellularTracks.erase(cellularTracks.begin()+removeTrack);
    }
  }
  
  return cellularTracks;
  
}

bool ConformalTracking::toBeUpdated(std::vector<cellularTrack> cellularTracks){
  bool update=false;
  for(unsigned int iTrack=0;iTrack<cellularTracks.size();iTrack++) if( cellularTracks[iTrack].back()->getWeight() > 0 ){update = true; break;}
  return update;
}

void ConformalTracking::followPath(std::vector<cellularTrack>& cellularTracks, int trackNumber, std::map<Cell*,bool>& usedCells, int& removeTrack){
  
  
  Cell* cell = cellularTracks[trackNumber].back();
  
  if(cell->getWeight() == 0) return;
  
  //  std::cout<<"  -- following path of cell with weight "<<cell->getWeight()<<std::endl;
  //  std::cout<<"  -- this cell is connected to "<<cell->getFrom().size()<<std::endl;
  if(cell->getFrom().size() > 1){
    //    std::cout<<"  -- branching!"<<std::endl;
    if(usedCells.count(cell->getFrom()[0])){
      removeTrack=trackNumber;
    }else{
      if( (cell->getWeight() - cell->getFrom()[0]->getWeight()) != 1  ){
        removeTrack=trackNumber;
      }else{
        cellularTracks[trackNumber].push_back(cell->getFrom()[0]);
        //        std::cout<<"  -- created branch 0"<<std::endl;
      }
    }
    for(unsigned int itCell=1;itCell<cell->getFrom().size();itCell++){
      if(usedCells.count(cell->getFrom()[itCell])){
        continue;
      }else{
        if( (cell->getWeight() - cell->getFrom()[itCell]->getWeight()) != 1  ) continue;
        //        std::cout<<"  -- created branch "<<itCell<<std::endl;
        //        std::cout<<"  -- cell weight is "<<cell->getWeight()<<" and branch has weight "<<cell->getFrom()[itCell]->getWeight()<<std::endl;
        cellularTrack branchedTrack = cellularTracks[trackNumber];
        branchedTrack.push_back(cell->getFrom()[itCell]);
        cellularTracks.push_back(branchedTrack);
      }
    }
  }
  
  while(cell->getWeight() > 0 && cell->getFrom().size() == 1){
    //    std::cout<<"  -- following linear path. Current weight "<<cell->getWeight()<<std::endl;
    Cell* parentCell = cell->getFrom()[0];
    if(usedCells.count(cell)){removeTrack=trackNumber; break;}
    cellularTracks[trackNumber].push_back(parentCell);
    cell = parentCell;
  }
  
  //	if(cell->getFrom().size() == 0) return;
  
  return;
  
}

// Given a list of connected cells (so-called cellular tracks), return the candidate(s) with lowest chi2/degrees of freedom.
// Attempt to remove hits to see if there is a significant improvement in the chi2/ndof, to protect against noise hits causing
// a good track to be discarded. If several candidates have low chi2/ndof and are not clones (limited sharing of hits) then
// return all of them. Given that there is no material scattering taken into account this helps retain low pt tracks, which
// may have worse chi2/ndof than ghosts/real tracks with an additional unrelated hit from the low pt track.
std::vector<KDTrack> ConformalTracking::getFittedTracks(std::vector<cellularTrack>& candidateTracks, std::vector<double>& finalChi2ndofs){
  
  // Make a container for all tracks being considered, initialise variables
  std::vector<KDTrack> trackContainer;
  std::vector<double> trackChi2ndofs;
  
  // Loop over all candidate tracks and do an inital fit to get the track angle (needed to calculate the
  // hit errors for the error-weighted fit)
  for(unsigned int itTrack=0;itTrack<candidateTracks.size();itTrack++){
    
		// Make the fitting object. TGraphErrors used for 2D error-weighted fitting
    KDTrack track;

    
    // Loop over all hits and add them to the fitter (and track)
    double npoints=0.;
    KDCluster* kdStart = candidateTracks[itTrack][0]->getEnd();
    track.add(kdStart); npoints++;
    
    for(unsigned int trackCell=0;trackCell<candidateTracks[itTrack].size();trackCell++){
      KDCluster* kdEnd = candidateTracks[itTrack][trackCell]->getStart();
      track.add(kdEnd); npoints++;
    }
    
    // Set up the track fitting
    ROOT::Math::Functor FCNFunction(track,2);
    newFitter.SetFunction(FCNFunction);
    newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
    newFitter.SetVariable(1,"intercept", 0., 0.1);
    newFitter.Minimize();
    
    // Fit the track
    track.fit(newFitter.X()[0],newFitter.X()[1]);
    double chi2ndof = track.chi2()/(npoints-2);
    
    // We try to see if there are spurious hits causing the chi2 to be very large. This would cause us to throw away
    // good tracks with perhaps just a single bad hit. Try to remove each hit and see if fitting without it causes a
    // significant improvement in chi2/ndof
    
    // Loop over each hit (starting at the back, since we will use the 'erase' function to get rid of them)
    // and see if removing it improves the chi2/ndof
    double removed=0;
    if(chi2ndof > 10. && chi2ndof < 300.){
      for(int point=npoints-1;point>=0;point--){
        
        // Stop if we would remove too many points on the track to meet the minimum hit requirement (or if the track has more
        // than 2 hits removed)
        if((npoints-removed-1) < m_minClustersOnTrack || removed == 2) break;
        
        // Refit the track without this point
        double newChi2ndof = fitWithoutPoint(track,point);
        
        // If the chi2/ndof is significantly better, remove the point permanently
        if( (chi2ndof - newChi2ndof) > 0 && (chi2ndof - newChi2ndof) > 0.5*chi2ndof ){
          track.remove(point);
          removed++;
          chi2ndof = newChi2ndof;
        }
      }
    }
  
    // Store the track information
    trackContainer.push_back(track);
    trackChi2ndofs.push_back(chi2ndof);
  }
  
  // Now have all sets of conformal tracks and their chi2/ndof. Decide which tracks to send back, ie. the one with
  // lowest chi2/ndof, and possibly others if they are not clones and have similar chi2 value
  std::vector<KDTrack> finalTracks = getLowestChi2(trackContainer,trackChi2ndofs,finalChi2ndofs);
  
  // Send back the final set of tracks
  return finalTracks;
  
}

// Pick the lowest chi2/ndof KDTrack from a list of possible tracks, and additionally return other tracks in the collection with similar chi2/ndof values that don't share many hits
std::vector<KDTrack> ConformalTracking::getLowestChi2(std::vector<KDTrack> trackContainer, std::vector<double> trackChi2ndofs, std::vector<double>& finalChi2ndofs){
  
  // Get the lowest chi2/ndof value from the given tracks
  double lowestChi2ndof = *std::min_element(trackChi2ndofs.begin(),trackChi2ndofs.end());
  KDTrack lowestChi2ndofTrack;
  for(unsigned int itTrack=0;itTrack<trackContainer.size();itTrack++) if(trackChi2ndofs[itTrack] == lowestChi2ndof) lowestChi2ndofTrack = trackContainer[itTrack];
  
	// Final track storage
  std::vector<KDTrack> finalTracks;
  
  // Loop over all other tracks and decide whether or not to save them
  for(unsigned int itTrack=0;itTrack<trackContainer.size();itTrack++){
    
    // Look at the difference in chi2/ndof - we want to keep tracks with similar chi2/ndof. If they
    // are clones then take the longest
    if( (trackChi2ndofs[itTrack] - lowestChi2ndof) < 10. ){
      
      // If same track and longer
      if(sameTrack(trackContainer[itTrack], lowestChi2ndofTrack)){
        if(trackContainer[itTrack].nPoints() > lowestChi2ndofTrack.nPoints()){
        	lowestChi2ndofTrack = trackContainer[itTrack];
        	lowestChi2ndof = trackChi2ndofs[itTrack];
        }
        continue;
      }
      
      // Store this track
      finalTracks.push_back(trackContainer[itTrack]);
      finalChi2ndofs.push_back(trackChi2ndofs[itTrack]);
    }
  }
  
  // Save the track with the lowest chi2/ndof
	finalTracks.insert(finalTracks.begin(),lowestChi2ndofTrack);
	finalChi2ndofs.insert(finalChi2ndofs.begin(),lowestChi2ndof);

  return finalTracks;
  
}

// Function to check if two KDtracks contain several hits in common
bool ConformalTracking::sameTrack(KDTrack track1, KDTrack track2){
  
  // Loop over all hits on track 1 and check if that hit is in track 2
  double nHitsInCommon = 0;
  std::vector<KDCluster*> track1hits = track1.clusters();
  std::vector<KDCluster*> track2hits = track2.clusters();
  
  for(int hit=0;hit<track1hits.size();hit++){
    if( std::find(track2hits.begin(),track2hits.end(),track1hits[hit]) !=  track2hits.end()) nHitsInCommon++;
  }
  
  // Cut on number of shared hits
  if(nHitsInCommon > 0.4 * track1hits.size()) return true;
  return false;
  
}

double ConformalTracking::fitWithoutPoint(KDTrack track,int point){

  // Remove the given point from the track
  track.remove(point);

  // Set up the fitter
  int npoints = track.nPoints();
  ROOT::Math::Functor FCNFunction(track,2);
  newFitter.SetFunction(FCNFunction);
  newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
  newFitter.SetVariable(1,"intercept", 0., 0.1);
  newFitter.Minimize();
  
  // Fit the track
  track.fit(newFitter.X()[0],newFitter.X()[1]);
  double chi2ndof = track.chi2ndof();
  
  return chi2ndof;
}

void ConformalTracking::updateCell(Cell* cell){
  
  if(cell->getTo().size() != 0){
    for(unsigned int i=0;i<cell->getTo().size();i++){
      cell->getTo()[i]->update(cell);
      updateCell(cell->getTo()[i]);
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

void ConformalTracking::extendTrack(KDTrack& track,std::vector<cellularTrack> trackSegments, std::map<KDCluster*,bool>& used){

  // Get the inital track chi2/ndof
  int npoints = track.nPoints();
  double chi2ndof = track.chi2ndof();

  // Of all of the track segments, get the one with lowest chi2/ndof. Now look at each point on the cellular
  // track and add it to the track. If the delta chi2/ndof is small enough, keep the hit
  std::vector<double> chi2ndofSegments;
  KDTrack bestTrackSegment = getFittedTracks(trackSegments,chi2ndofSegments)[0];
  
  double newChi2ndof = fitWithExtension(track, bestTrackSegment.clusters());

//  if(fabs(newChi2ndof-chi2ndof) < 2.*chi2ndof || (newChi2ndof-chi2ndof) < 10.){
  if( newChi2ndof < m_chi2cut ){
    for(int newpoint=(bestTrackSegment.clusters().size()-3);newpoint>=0;newpoint--){
      track.insert(bestTrackSegment.clusters()[newpoint]);
      if(newChi2ndof < 10.) used[bestTrackSegment.clusters()[newpoint]] = true;
      chi2ndof = newChi2ndof;
    }
  }
  
}

double ConformalTracking::fitWithExtension(KDTrack track, std::vector<KDCluster*> hits){
  
  // Add the point to the track
  for(int i=(hits.size()-3);i>=0;i--) track.add(hits[i]);

  int npoints = track.nPoints();
  
  // Fit the track and get the new chi2
  ROOT::Math::Functor FCNFunction(track,2);
  newFitter.SetFunction(FCNFunction);
  newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
  newFitter.SetVariable(1,"intercept", 0., 0.1);
  newFitter.Minimize();
  
  track.fit(newFitter.X()[0],newFitter.X()[1]);
  double chi2ndof = track.chi2ndof();
  
  return chi2ndof;

}

double ConformalTracking::fitWithPoint(KDTrack track, KDCluster* hit){

  // Add the point to the track
  track.add(hit);
  int npoints = track.nPoints();
  
  // Fit the track and get the new chi2
  ROOT::Math::Functor FCNFunction(track,2);
  newFitter.SetFunction(FCNFunction);
  newFitter.SetVariable(0,"gradient", track.clusters()[npoints-1]->getV()/track.clusters()[npoints-1]->getU(), 0.1);
  newFitter.SetVariable(1,"intercept", 0., 0.1);
  newFitter.Minimize();
  
  track.fit(newFitter.X()[0],newFitter.X()[1]);
  double chi2ndof = track.chi2ndof();

  return chi2ndof;

}


// Debug function - checks if a track will be associated to an MC particle or not
double ConformalTracking::checkReal(cellularTrack track, std::map<KDCluster*,MCParticle*> kdParticles, std::map<MCParticle*,bool>& reconstructed){
 
  // Store all mcparticles associated to this track
  std::vector<MCParticle*> particles;
  std::map<MCParticle*,double> particleHits;
  double nHits=0.;
  
 	// Get the first hit on the track
  KDCluster* kdEnd = track[0]->getEnd();
  if(!kdEnd->removed()){
  	nHits++;
  
    // Fill the MC particle info
    if(particleHits.count(kdParticles[kdEnd])){
      particleHits[kdParticles[kdEnd]]++;
    }else{
      particles.push_back(kdParticles[kdEnd]);
      particleHits[kdParticles[kdEnd]] = 1;
    }
    
  }

  // Loop over all cells and get the hit that they connect to
  for(unsigned int trackCell=0;trackCell<track.size();trackCell++){
    KDCluster* kdStart = track[trackCell]->getStart();
    if(kdStart->removed())continue;
    nHits++;
    // Fill the MC particle info
    if(particleHits.count(kdParticles[kdStart])){
      particleHits[kdParticles[kdStart]]++;
    }else{
      particles.push_back(kdParticles[kdStart]);
      particleHits[kdParticles[kdStart]] = 1.;
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
  std::cout<<"Number of hits on track: "<<nHits<<". Good hits: "<<bestHits<<". Purity: "<<purity<<". Pt: "<<sqrt( bestParticle->getMomentum()[0]*bestParticle->getMomentum()[0] + bestParticle->getMomentum()[1]*bestParticle->getMomentum()[1] )<<std::endl;
  
  if(purity >= 0.75){
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

      if(hit->sameLayer(nhit)) sameID = false;
    }
    
    // Didn't find another hit with the same ID
    if(!sameID) nUniqueHits++;
    
  }
  
  return nUniqueHits;
}






