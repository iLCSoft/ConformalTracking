#ifndef PatternRecognition_h
#define PatternRecognition_h 1

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
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"

#include <AIDA/AIDA.h>

#include "KDTree.h"
#include "Cell.h"

using namespace lcio ;
using namespace marlin ;
using namespace AIDA ;

class PatternRecognition : public Processor {
		
public:
	
	virtual Processor*  newProcessor() { return new PatternRecognition ; }
	
	PatternRecognition() ;
	
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
	void drawline(KDCluster*, KDCluster*, int);
	
	// Pattern recognition algorithms
	std::vector<cellularTrack> createTracks(Cell*, std::map<Cell*, bool>&);
	bool toBeUpdated(std::vector<cellularTrack>);
	void followPath(std::vector<cellularTrack>&, int, std::map<Cell*, bool>&, int&);
  cellularTrack getLowestChi2(std::vector<cellularTrack>);
  void updateCell(Cell*);
  KDCluster* extrapolateCell(Cell*, double);
  
  void extendSeedCells(std::vector<Cell*>&, std::map<KDCluster*,bool>, KDTree*, bool);

  
	
protected:
	
	// Collection names for (in/out)put
	std::vector<std::string> m_inputTrackerHitCollections ;
	std::string m_outputTrackCollection ;
	
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
	TH1F* m_theta;
	TH1F* m_cellAngle;
	TH2F* m_cellAngleRadius;
	TH2F* m_cellAngleInverseRadius;
	TH2F* m_conformalEvents;
	TH2F* m_conformalEventsZU;
	TH2F* m_conformalEventsZV;
	TH2F* m_nonconformalEvents;
	TH2F* m_conformalEventsRTheta;
	TH3F* m_conformalEvents3D;
	
	TH1F* m_cellAngleMC;
	TH2F* m_cellAngleRadiusMC;
	TH2F* m_cellAngleInverseRadiusMC;
	TH2F* m_gradientRadiusMC;
	
	// Other constants
  double m_thetaRange;
  double m_maxCellAngle;
  double m_maxDistance;
  int m_minClustersOnTrack;
  bool m_debugPlots;
	
} ;

bool sort_by_radius(EVENT::TrackerHit*, EVENT::TrackerHit*);
bool sort_by_radiusKD(KDCluster*, KDCluster*);
bool sort_by_cellWeight(Cell*, Cell*);

#endif



