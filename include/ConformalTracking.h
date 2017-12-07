#ifndef ConformalTracking_h
#define ConformalTracking_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>
//#include "DDRec/Surface.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include "EVENT/TrackerHit.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include <AIDA/AIDA.h>

#include "Cell.h"
#include "DebugTool.h"
#include "KDTrack.h"
#include "KDTree.h"
#include "KalmanTrack.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

using namespace lcio;
using namespace marlin;
using namespace AIDA;

class ConformalTracking : public Processor {
  struct Parameters {
  public:
    double _maxCellAngle;
    double _maxCellAngleRZ;
    double _chi2cut;
    int    _minClustersOnTrack;
    double _maxDistance;
    bool   _highPTfit;
    bool   _onlyZSchi2cut;
    Parameters(double maxCellAngle, double maxCellAngleRZ, double chi2cut, int minClustersOnTrack, double maxDistance,
               bool highPTfit, bool onlyZSchi2cut)
        : _maxCellAngle(maxCellAngle),
          _maxCellAngleRZ(maxCellAngleRZ),
          _chi2cut(chi2cut),
          _minClustersOnTrack(minClustersOnTrack),
          _maxDistance(maxDistance),
          _highPTfit(highPTfit),
          _onlyZSchi2cut(onlyZSchi2cut) {}
    Parameters(Parameters const& rhs) = default;
  };

public:
  virtual Processor* newProcessor() { return new ConformalTracking; }

  ConformalTracking();
  ConformalTracking(const ConformalTracking&) = delete;
  ConformalTracking& operator=(const ConformalTracking&) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init();

  // Called at the beginning of every run
  virtual void processRunHeader(LCRunHeader*) { m_runNumber++; }

  // Run over each event - the main algorithm
  virtual void processEvent(LCEvent* evt);

  // Run at the end of each event
  virtual void check(LCEvent*){};

  // Called at the very end for cleanup, histogram saving, etc.
  virtual void end();

  // Call to get collections
  void getCollection(LCCollection*&, std::string const&, LCEvent*);

  // Plotting function for displaying cells
  void drawline(SKDCluster const&, SKDCluster const&, int, int style = 1);

  // Pattern recognition algorithms:

  // Cell creation
  SKDCluster extrapolateCell(Cell::SCell const&, double);
  void       extendSeedCells(SharedCells&, UKDTree&, bool, const SharedKDClusters&, Parameters const&,
                       bool vertexToTracker = true);

  // Track finding
  void buildNewTracks(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool radialSearch = false,
                      bool vertexToTracker = true);
  void extendTracks(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&);
  void combineCollections(SharedKDClusters&, UKDTree&, std::vector<int> const&, std::map<int, SharedKDClusters> const&);

  void extendHighPT(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool radialSearch = false);

  void createTracksNew(UniqueCellularTracks&, Cell::SCell&, std::map<Cell::SCell, bool>&);
  bool toBeUpdated(UniqueCellularTracks const&);
  void updateCell(Cell::SCell const&);

  // Track fitting
  void getFittedTracks(UniqueKDTracks&, UniqueCellularTracks&, std::map<Cell::SCell, bool>&, Parameters const&);
  void getLowestChi2(UniqueKDTracks&, UniqueKDTracks&);

  double fitWithoutPoint(KDTrack, int);
  int    overlappingHits(const UKDTrack&, const UKDTrack&);

  void extendTrack(UKDTrack&, UniqueCellularTracks, std::map<SKDCluster, bool>&, std::map<Cell::SCell, bool>&);
  //double fitWithPoint(KalmanTrack, KDCluster*);
  void fitWithPoint(KDTrack, SKDCluster&, double&, double&);

  double fitWithExtension(KDTrack, SharedKDClusters, double&, double&);

  // MC truth debug
  double checkReal(UKDTrack&, std::map<SKDCluster, MCParticle*>, std::map<MCParticle*, bool>&,
                   std::map<MCParticle*, SharedKDClusters>);
  int  getUniqueHits(SharedKDClusters);
  void checkReconstructionFailure(MCParticle*, std::map<MCParticle*, SharedKDClusters>, UKDTree&, Parameters const&);
  void checkUnallowedTracks(UniqueCellularTracks, Parameters const&);

protected:
  // Collection names for (in/out)put
  std::vector<std::string> m_inputTrackerHitCollections{};
  std::string              m_outputTrackCollection{};
  std::string              m_inputParticleCollection{};
  std::vector<std::string> m_inputRelationCollections{};
  std::string              m_outputDebugHits{};

  // Run and event counters
  int m_eventNumber = 0;
  int m_runNumber   = 0;

  // Track fit factory
  MarlinTrk::IMarlinTrkSystem* trackFactory = nullptr;

  // Track fit parameters
  double m_initialTrackError_d0    = 0.0;
  double m_initialTrackError_phi0  = 0.0;
  double m_initialTrackError_omega = 0.0;
  double m_initialTrackError_z0    = 0.0;
  double m_initialTrackError_tanL  = 0.0;
  double m_maxChi2perHit           = 0.0;
  double m_magneticField           = 0.0;

  // Histograms
  TH1F* m_cellAngle           = nullptr;
  TH1F* m_cellDOCA            = nullptr;
  TH2F* m_cellAngleRadius     = nullptr;
  TH2F* m_cellLengthRadius    = nullptr;
  TH2F* m_cellAngleLength     = nullptr;
  TH1F* m_conformalChi2       = nullptr;
  TH1F* m_conformalChi2real   = nullptr;
  TH1F* m_conformalChi2fake   = nullptr;
  TH2F* m_conformalChi2Purity = nullptr;

  TH1F* m_absZ = nullptr;

  TH1F* m_conformalChi2MC        = nullptr;
  TH2F* m_conformalChi2PtMC      = nullptr;
  TH2F* m_conformalChi2VertexRMC = nullptr;

  TH1F* m_conformalChi2SzMC        = nullptr;
  TH2F* m_conformalChi2SzPtMC      = nullptr;
  TH2F* m_conformalChi2SzVertexRMC = nullptr;

  TH1F* m_cellAngleMC        = nullptr;
  TH1F* m_cellDOCAMC         = nullptr;
  TH1F* m_cellAngleRZMC      = nullptr;
  TH2F* m_cellAngleRadiusMC  = nullptr;
  TH2F* m_cellLengthRadiusMC = nullptr;
  TH2F* m_cellAngleLengthMC  = nullptr;

  TH2F* m_conformalEvents       = nullptr;
  TH2F* m_nonconformalEvents    = nullptr;
  TH2F* m_conformalEventsRTheta = nullptr;
  TH2F* m_conformalEventsMC     = nullptr;

  TCanvas* m_canvConformalEventDisplay                  = nullptr;
  TCanvas* m_canvConformalEventDisplayAllCells          = nullptr;
  TCanvas* m_canvConformalEventDisplayAcceptedCells     = nullptr;
  TCanvas* m_canvConformalEventDisplayMC                = nullptr;
  TCanvas* m_canvConformalEventDisplayMCunreconstructed = nullptr;

  TH2F* m_szDistribution  = nullptr;
  TH2F* m_uvDistribution  = nullptr;
  TH2F* m_xyDistribution  = nullptr;
  TH3F* m_xyzDistribution = nullptr;

  // Other constants
  double            m_thetaRange         = 0.0;
  double            m_chi2cut            = 0.0;
  double            m_chi2increase       = 0.0;
  double            m_maxCellAngle       = 0.0;
  double            m_maxCellAngleRZ     = 0.0;
  double            m_maxDistance        = 0.0;
  int               m_minClustersOnTrack = 0;
  bool              m_debugPlots         = false;
  bool              m_retryTooManyTracks = true;
  bool              m_sortTreeResults    = true;
  double            m_purity             = 0.0;
  SKDCluster        debugSeed            = nullptr;
  ConformalDebugger m_debugger{};
};

// ---------------------------
// SORT FUNCTIONS
// ---------------------------

// Sort tracker hits from smaller to larger radius
bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2) {
  double radius1 =
      sqrt((hit1->getPosition()[0]) * (hit1->getPosition()[0]) + (hit1->getPosition()[1]) * (hit1->getPosition()[1]));
  double radius2 =
      sqrt((hit2->getPosition()[0]) * (hit2->getPosition()[0]) + (hit2->getPosition()[1]) * (hit2->getPosition()[1]));
  return (radius1 < radius2);
}

// Sort kd hits from larger to smaller radius
bool sort_by_radiusKD(SKDCluster const& hit1, SKDCluster const& hit2) {
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 > radius2);
}

// Sort kd hits from smaller to larger radius
bool sort_by_lower_radiusKD(SKDCluster const& hit1, SKDCluster const& hit2) {
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 < radius2);
}

// Sort kdhits by lower to higher layer number
bool sort_by_layer(SKDCluster const& hit1, SKDCluster const& hit2) {
  if (hit1->getSubdetector() != hit2->getSubdetector())
    return (hit1->getSubdetector() < hit2->getSubdetector());
  else if (hit1->getLayer() != hit2->getLayer())
    return (hit1->getLayer() < hit2->getLayer());
  else if (hit1->getSide() != hit2->getSide())
    return (hit1->getSide() < hit2->getSide());
  else
    return false;
}

// Sort cells from higher to lower weight
bool sort_by_cellWeight(Cell::SCell const& cell1, Cell::SCell const& cell2) {
  int weight1 = cell1->getWeight();
  int weight2 = cell2->getWeight();
  return (weight1 > weight2);
}

// Sort kdtracks from longest to shortest
bool sort_by_length(UKDTrack const& track1, UKDTrack const& track2) {
  return (track1->m_clusters.size() > track2->m_clusters.size());
}

// Sort kdtracks from lowest to highest pt
bool sort_by_pt(UKDTrack const& track1, UKDTrack const& track2) { return (track1->pt() > track2->pt()); }
#endif
