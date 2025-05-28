#ifndef ConformalTracking_h
#define ConformalTracking_h 1

#include "marlin/Processor.h"

#include "Cell.h"
#include "DebugTool.h"
#include "KDTrack.h"
#include "KDTree.h"
#include "Parameters.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include "lcio.h"

#include <MarlinTrk/IMarlinTrkSystem.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using namespace marlin;

class ConformalTracking : public Processor {
public:
  virtual Processor* newProcessor() { return new ConformalTracking; }

  ConformalTracking();
  ConformalTracking(std::string const& procName);
  ConformalTracking(const ConformalTracking&)            = delete;
  ConformalTracking& operator=(const ConformalTracking&) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init();

  // Register input collections and parameters
  void registerParameters();

  /// fill the vectors mapping collections to index
  void fillCollectionIndexVectors();

  // Called at the beginning of every run
  virtual void processRunHeader(LCRunHeader*) { m_runNumber++; }

  // Run over each event - the main algorithm
  virtual void processEvent(LCEvent* evt);

  // Run at the end of each event
  virtual void check(LCEvent*) {};

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

  void runStep(SharedKDClusters&, UKDTree&, UniqueKDTracks&, std::map<int, SharedKDClusters> const&, Parameters const&);
  virtual void parseStepParameters();

  // Track finding
  void buildNewTracks(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool radialSearch = false,
                      bool vertexToTracker = true);
  bool neighbourIsCompatible(const SKDCluster& neighbourHit, const SKDCluster& seedHit, const double slopeZRange);
  void extendTracks(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&);
  void combineCollections(SharedKDClusters&, UKDTree&, std::vector<int> const&, std::map<int, SharedKDClusters> const&);

  void extendHighPT(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool radialSearch = false);

  void extendTracksPerLayer(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool vertexToTracker = true);

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
  double checkReal(UKDTrack&, std::map<MCParticle*, bool>&, std::map<MCParticle*, SharedKDClusters>);
  int    getUniqueHits(SharedKDClusters);
  void   checkReconstructionFailure(MCParticle*, std::map<MCParticle*, SharedKDClusters>, UKDTree&, Parameters const&);
  void   checkUnallowedTracks(UniqueCellularTracks, Parameters const&);

protected:
  std::vector<Parameters> _stepParameters{};

  // Collection names for (in/out)put
  std::vector<std::string> m_inputTrackerHitCollections{};
  std::vector<std::string> m_inputMainTrackerHitCollections{};
  std::vector<std::string> m_inputVertexBarrelCollections{};
  std::vector<std::string> m_inputVertexEndcapCollections{};

  std::string              m_outputTrackCollection{};
  std::string              m_inputParticleCollection{};
  std::vector<std::string> m_inputRelationCollections{};
  std::string              m_outputDebugHits{};
  std::vector<int>         m_allHits{};
  std::vector<int>         m_trackerHits{};
  std::vector<int>         m_vertexBarrelHits{};
  std::vector<int>         m_vertexEndcapHits{};
  std::vector<int>         m_vertexCombinedHits{};

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
  TH1F* m_X                 = nullptr;
  TH1F* m_Y                 = nullptr;
  TH1F* m_Z                 = nullptr;
  TH1F* m_neighX            = nullptr;
  TH1F* m_neighY            = nullptr;
  TH1F* m_neighZ            = nullptr;
  TH1F* m_slopeZ            = nullptr;
  TH1F* m_slopeZ_true       = nullptr;
  TH1F* m_slopeZ_true_first = nullptr;
  TH2F* m_slopeZ_vs_pt_true = nullptr;

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
  double                               m_thetaRange                 = 0.0;
  double                               m_chi2cut                    = 0.0;
  double                               m_maxCellAngle               = 0.0;
  double                               m_maxCellAngleRZ             = 0.0;
  double                               m_maxDistance                = 0.0;
  double                               m_slopeZRange                = 1000.0;
  double                               m_highPTcut                  = 0.0;
  int                                  m_minClustersOnTrack         = 0;
  int                                  m_minClustersOnTrackAfterFit = 0;
  int                                  m_maxHitsInvFit              = 0;
  int                                  m_tooManyTracks              = 5e5;
  bool                                 m_enableTCVC                 = true;
  bool                                 m_debugPlots                 = false;
  bool                                 m_debugTime                  = false;
  bool                                 m_retryTooManyTracks         = true;
  bool                                 m_sortTreeResults            = true;
  double                               m_purity                     = 0.0;
  SKDCluster                           debugSeed                    = nullptr;
  ConformalDebugger                    m_debugger{};
  std::map<SKDCluster, MCParticle*>    kdParticles{};  // Link from conformal hit to MC particle
  std::map<SKDCluster, SimTrackerHit*> kdSimHits{};    // Link from conformal hit to SimHit
};

// ---------------------------
// SORT FUNCTIONS
// ---------------------------

// Sort tracker hits from smaller to larger radius
inline bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2) {
  double radius1 =
      sqrt((hit1->getPosition()[0]) * (hit1->getPosition()[0]) + (hit1->getPosition()[1]) * (hit1->getPosition()[1]));
  double radius2 =
      sqrt((hit2->getPosition()[0]) * (hit2->getPosition()[0]) + (hit2->getPosition()[1]) * (hit2->getPosition()[1]));
  return (radius1 < radius2);
}

// Sort kd hits from larger to smaller radius
inline bool sort_by_radiusKD(SKDCluster const& hit1, SKDCluster const& hit2) {
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 > radius2);
}

// Sort kd hits from smaller to larger radius
inline bool sort_by_lower_radiusKD(SKDCluster const& hit1, SKDCluster const& hit2) {
  double radius1 = hit1->getR();
  double radius2 = hit2->getR();
  return (radius1 < radius2);
}

// Sort kdhits by lower to higher layer number
inline bool sort_by_layer(SKDCluster const& hit1, SKDCluster const& hit2) {
  if (hit1->getSubdetector() != hit2->getSubdetector())
    return (hit1->getSubdetector() < hit2->getSubdetector());
  else if (hit1->getLayer() != hit2->getLayer())
    return (hit1->getLayer() < hit2->getLayer());
  else if (hit1->getSide() != hit2->getSide())
    return (hit1->getSide() < hit2->getSide());
  else
    return false;
}

// Sort kdhits by higher to lower layer number
inline bool sort_by_lower_layer(SKDCluster const& hit1, SKDCluster const& hit2) {
  if (hit1->getSubdetector() != hit2->getSubdetector())
    return (hit1->getSubdetector() > hit2->getSubdetector());
  else if (hit1->getLayer() != hit2->getLayer())
    return (hit1->getLayer() > hit2->getLayer());
  else if (hit1->getSide() != hit2->getSide())
    return (hit1->getSide() > hit2->getSide());
  else
    return false;
}

// Sort cells from higher to lower weight
inline bool sort_by_cellWeight(Cell::SCell const& cell1, Cell::SCell const& cell2) {
  int weight1 = cell1->getWeight();
  int weight2 = cell2->getWeight();
  return (weight1 > weight2);
}

// Sort kdtracks from longest to shortest
inline bool sort_by_length(UKDTrack const& track1, UKDTrack const& track2) {
  return (track1->m_clusters.size() > track2->m_clusters.size());
}

// Sort kdtracks from lowest to highest pt
inline bool sort_by_pt(UKDTrack const& track1, UKDTrack const& track2) { return (track1->pt() > track2->pt()); }
#endif
