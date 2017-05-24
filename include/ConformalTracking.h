#ifndef ConformalTracking_h
#define ConformalTracking_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <EVENT/LCCollection.h>

#include <gsl/gsl_rng.h>
//#include "DDRec/Surface.h"
#include <EVENT/LCCollection.h>
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "EVENT/TrackerHit.h"
#include <EVENT/MCParticle.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

#include <AIDA/AIDA.h>

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "KDTrack.h"
#include "KDTree.h"
#include "Cell.h"
#include "KalmanTrack.h"

using namespace lcio ;
using namespace marlin ;
using namespace AIDA ;

class ConformalTracking : public Processor {
    
public:
    
    virtual Processor*  newProcessor() { return new ConformalTracking ; }
    
    ConformalTracking() ;
    
    // Initialisation - run at the beginning to start histograms, etc.
    virtual void init() ;
    
    // Called at the beginning of every run
    virtual void processRunHeader(LCRunHeader* run){m_runNumber++;}
    
    // Run over each event - the main algorithm
    virtual void processEvent(LCEvent* evt) ;
    
    // Run at the end of each event
    virtual void check(LCEvent* evt){} ;
    
    // Called at the very end for cleanup, histogram saving, etc.
    virtual void end() ;
    
    
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);
    
    // Plotting function for displaying cells
    void drawline(KDCluster*, KDCluster*, int, int style=1);
    
    // Pattern recognition algorithms:
    
    // Cell creation
    KDCluster* extrapolateCell(Cell*, double);
    void extendSeedCells(std::vector<Cell*>&, KDTree*, bool, std::vector<KDCluster*>);
    
    // Track finding
  void buildNewTracks(std::vector<KDTrack*>&, std::vector<KDCluster*>&, KDTree*);
  void extendTracks(std::vector<KDTrack*>&, std::vector<KDCluster*>&, KDTree*);
  void combineCollections(std::vector<KDCluster*>&, KDTree*&, std::vector<int>, std::map<int,std::vector<KDCluster*>>);
  
    void createTracksNew(std::vector<cellularTrack*>&, Cell*, std::map<Cell*, bool>&);
    bool toBeUpdated(std::vector<cellularTrack*>const&);
    void updateCell(Cell*);
    
    // Track fitting
    void getFittedTracks(std::vector<KDTrack*>&, std::vector<cellularTrack*>&, std::map<Cell*,bool>&);
    void getLowestChi2(std::vector<KDTrack*>&, std::vector<KDTrack*>);
    
    double fitWithoutPoint(KDTrack,int);
    int overlappingHits(const KDTrack*, const KDTrack*);
    
    void extendTrack(KDTrack*,std::vector<cellularTrack*>,std::map<KDCluster*,bool>&, std::map<Cell*,bool>&);
    double fitWithPoint(KalmanTrack, KDCluster*);
    double fitWithPoint(KDTrack, KDCluster*);
    
    double fitWithExtension(KDTrack, std::vector<KDCluster*>, double&, double&);
    
    // MC truth debug
    double checkReal(KDTrack*, std::map<KDCluster*,MCParticle*>, std::map<MCParticle*,bool>&, std::map<MCParticle*, std::vector<KDCluster*> >);
    int getUniqueHits(std::vector<KDCluster*>);
    void checkReconstructionFailure(MCParticle*, std::map<MCParticle*, std::vector<KDCluster*> >, KDTree*);
    
    
protected:
    
    // Collection names for (in/out)put
    std::vector<std::string> m_inputTrackerHitCollections ;
    std::string m_outputTrackCollection ;
    std::string m_inputParticleCollection;
    std::vector<std::string> m_inputRelationCollections ;
    std::string m_outputDebugHits;
    
    // Run and event counters
    int m_eventNumber ;
    int m_runNumber ;
    
    // Track fit factory
    MarlinTrk::IMarlinTrkSystem* trackFactory;
    
    // Track fit parameters
    double m_initialTrackError_d0;
    double  m_initialTrackError_phi0;
    double m_initialTrackError_omega;
    double m_initialTrackError_z0;
    double m_initialTrackError_tanL;
    double m_maxChi2perHit;
    double m_magneticField;
    
    // Histograms
    TH1F* m_cellAngle;
    TH2F* m_cellAngleRadius;
    TH2F* m_cellLengthRadius;
    TH2F* m_cellAngleLength;
    TH1F* m_conformalChi2;
    TH1F* m_conformalChi2real;
    TH1F* m_conformalChi2fake;
    TH2F* m_conformalChi2Purity;
  
  TH1F* m_conformalChi2MC;
  TH2F* m_conformalChi2PtMC;
  TH2F* m_conformalChi2VertexRMC;

  TH1F* m_conformalChi2SzMC;
  TH2F* m_conformalChi2SzPtMC;
  TH2F* m_conformalChi2SzVertexRMC;

    TH1F* m_cellAngleMC;
    TH1F* m_cellAngleRZMC;
    TH2F* m_cellAngleRadiusMC;
    TH2F* m_cellLengthRadiusMC;
    TH2F* m_cellAngleLengthMC;
    
    TH2F* m_conformalEvents;
    TH2F* m_nonconformalEvents;
    TH2F* m_conformalEventsRTheta;
    TH2F* m_conformalEventsMC;
    
    TCanvas* m_canvConformalEventDisplay;
    TCanvas* m_canvConformalEventDisplayAllCells;
    TCanvas* m_canvConformalEventDisplayAcceptedCells;
    TCanvas* m_canvConformalEventDisplayMC;
    TCanvas* m_canvConformalEventDisplayMCunreconstructed;
    
    TH2F* m_szDistribution;
    TH2F* m_uvDistribution;
    TH2F* m_xyDistribution;
    TH3F* m_xyzDistribution;
    
    // Other constants
    double m_thetaRange;
    double m_chi2cut;
    double m_maxCellAngle;
    double m_maxCellAngleRZ;
    double m_maxDistance;
    int m_minClustersOnTrack;
    bool m_debugPlots;
    double m_purity;
    KDCluster* debugSeed;
    
} ;

// ---------------------------
// SORT FUNCTIONS
// ---------------------------

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

#endif



