#ifndef ConformalTracking_h
#define ConformalTracking_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <map>
#include <string>
#include <vector>

#include <EVENT/LCCollection.h>

#include <EVENT/LCCollection.h>
#include <gsl/gsl_rng.h>
#include "DDRec/Surface.h"
#include "EVENT/TrackerHit.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include <AIDA/AIDA.h>

#include "Cell.h"
#include "KDTree.h"

using namespace lcio;
using namespace marlin;
using namespace AIDA;

class ConformalTracking : public Processor {
public:
  virtual Processor* newProcessor() { return new ConformalTracking; }

  ConformalTracking();

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init();

  // Called at the beginning of every run
  virtual void processRunHeader(LCRunHeader* run);

  // Run over each event - the main algorithm
  virtual void processEvent(LCEvent* evt);

  // Run at the end of each event
  virtual void check(LCEvent* evt);

  // Called at the very end for cleanup, histogram saving, etc.
  virtual void end();

  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);

  // Plotting function for displaying cells
  void drawline(KDCluster*, KDCluster*, int);

  // Pattern recognition algorithms
  std::vector<cellularTrack> createTracks(Cell*, std::map<Cell*, bool>&);
  bool                       toBeUpdated(std::vector<cellularTrack>);
  void                       followPath(std::vector<cellularTrack>&, int, std::map<Cell*, bool>&, int&);
  cellularTrack              getLowestChi2(std::vector<cellularTrack>);
  void                       updateCell(Cell*);
  KDCluster*                 extrapolateCell(Cell*, double);

  void extendSeedCells(std::vector<Cell*>&, std::map<KDCluster*, bool>, KDTree*, bool);

protected:
  // Collection names for (in/out)put
  std::vector<std::string> m_inputTrackerHitCollections;
  std::string              m_outputTrackCollection;

  // Run and event counters
  int m_eventNumber;
  int m_runNumber;

  // Track fit factory
  MarlinTrk::IMarlinTrkSystem* trackFactory;

  // Track fit parameters
  double m_initialTrackError_d0;
  double m_initialTrackError_phi0;
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

  TH1F* m_cellAngleMC;
  TH2F* m_cellAngleRadiusMC;
  TH2F* m_cellLengthRadiusMC;
  TH2F* m_cellAngleLengthMC;

  TH2F* m_conformalEvents;
  TH2F* m_nonconformalEvents;
  TH2F* m_conformalEventsRTheta;
  TH2F* m_conformalEventsMC;

  TCanvas* m_canvConformalEventDisplay;
  TCanvas* m_canvConformalEventDisplayMC;

  // Other constants
  double m_thetaRange;
  double m_maxCellAngle;
  double m_maxDistance;
  int    m_minClustersOnTrack;
  bool   m_debugPlots;
};

bool sort_by_radius(EVENT::TrackerHit*, EVENT::TrackerHit*);
bool sort_by_radiusKD(KDCluster*, KDCluster*);
bool sort_by_cellWeight(Cell*, Cell*);

#endif
