#ifndef ConformalTracking_h
#define ConformalTracking_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include <EVENT/LCCollection.h>

#include <gsl/gsl_rng.h>
#include "DDRec/Surface.h"
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
	virtual void processRunHeader( LCRunHeader* run ) ;
	
	// Run over each event - the main algorithm
	virtual void processEvent( LCEvent * evt ) ;
	
  // Run at the end of each event
	virtual void check( LCEvent * evt ) ;
	
	// Called at the very end for cleanup, histogram saving, etc.
	virtual void end() ;
	
  
  
  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);

	// Plotting function for displaying cells
  void drawline(KDCluster*, KDCluster*, int, int style=1);
	
	// Pattern recognition algorithms:
  
  // Cell creation
  KDCluster* extrapolateCell(Cell*, double);
  void extendSeedCells(std::vector<Cell*>&, std::map<KDCluster*,bool>, KDTree*, bool, std::vector<KDCluster*>);

  // Track finding
  std::vector<cellularTrack> createTracksNew(Cell*, std::map<Cell*, bool>&);
	bool toBeUpdated(std::vector<cellularTrack>);
  void updateCell(Cell*);
  
  // Track fitting
  std::vector<KDTrack> getFittedTracks(std::vector<cellularTrack>&, std::vector<double>&, std::map<Cell*,bool>&);
  std::vector<KDTrack> getLowestChi2(std::vector<KDTrack>, std::vector<double>, std::vector<double>&);

  double chi2SZ(KDTrack);
  
  double fitWithoutPoint(KDTrack,int);
  int overlappingHits(KDTrack, KDTrack);
  
  void extendTrack(KDTrack&,std::vector<cellularTrack>,std::map<KDCluster*,bool>&, std::map<Cell*,bool>&);
  double fitWithPoint(KDTrack, KDCluster*);

  double fitWithExtension(KDTrack, std::vector<KDCluster*>);

  ROOT::Minuit2::Minuit2Minimizer newFitter;
//  ROOT::Math::Functor* FCNFunction;

  // MC truth debug
  double checkReal(KDTrack, std::map<KDCluster*,MCParticle*>, std::map<MCParticle*,bool>&);
  int getUniqueHits(std::vector<KDCluster*>);
  
	
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

  TH1F* m_cellAngleMC;
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
  double m_maxDistance;
  int m_minClustersOnTrack;
  bool m_debugPlots;
  KDCluster* debugSeed;
	
} ;

bool sort_by_radius(EVENT::TrackerHit*, EVENT::TrackerHit*);
bool sort_by_radiusKD(KDCluster*, KDCluster*);
bool sort_by_cellWeight(Cell*, Cell*);

#endif



