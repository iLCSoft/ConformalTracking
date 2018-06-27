#ifndef ConformalTrackingV2_h
#define ConformalTrackingV2_h 1

#include "ConformalTracking.h"

class ConformalTrackingV2 : public ConformalTracking {
public:
  virtual Processor* newProcessor() { return new ConformalTrackingV2; }

  ConformalTrackingV2();
  ConformalTrackingV2(const ConformalTrackingV2&) = delete;
  ConformalTrackingV2& operator=(const ConformalTrackingV2&) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init();

  // Called at the beginning of every run
  virtual void processRunHeader(LCRunHeader*) { m_runNumber++; }

  virtual void parseStepParameters();

protected:
  // clang-format off
  std::vector<std::string> m_rawSteps = {R"RAW(
    [VXDBarrel]
    @Collections = VXDTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, BuildNewTracks
    [VXDEncap]
    @Collections = VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, ExtendTracks
    [TCVC]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, BuildNewTracks
    [LowerCellAngle1]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.175; MaxCellAngleRZ = 0.175;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = CombineCollections, BuildNewTracks
    [LowerCellAngle2]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 6000;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = BuildNewTracks
    [LowerHitNumbers]
    @Collections = VXDTrackerHits,VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 6000;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = BuildNewTracks, SortTracks
    [Tracker]
    @Collections = ITrackerHits,OTrackerHits,ITrackerEndcapHits,OTrackerEndcapHits
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 6000;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = CombineCollections, ExtendTracks
    [Displaced]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits, ITrackerHits, OTrackerHits, ITrackerEndcapHits, OTrackerEndcapHits
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 3000;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = OnlyZSchi2cut, VertexToTracker, RadialSearch
    @Functions = CombineCollections, ExtendTracks
)RAW"
  };
  // clang-format on
};

#endif  // ConformalTrackingV2_h
