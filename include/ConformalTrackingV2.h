#ifndef ConformalTrackingV2_h
#define ConformalTrackingV2_h 1

#include "ConformalTracking.h"

class ConformalTrackingV2 : public ConformalTracking {
public:
  virtual Processor* newProcessor() { return new ConformalTrackingV2; }

  ConformalTrackingV2();
  ConformalTrackingV2(const ConformalTrackingV2&) = delete;
  ConformalTrackingV2& operator=(const ConformalTrackingV2&) = delete;

  // Called at the beginning of every run
  virtual void processRunHeader(LCRunHeader*) { m_runNumber++; }

  virtual void parseStepParameters();

protected:
  // clang-format off
  std::vector<std::string> m_rawSteps = {R"RAW(
    [VXDBarrel]
    @Collections : VXDTrackerHits
    @Parameters : MaxCellAngle : 0.005; MaxCellAngleRZ : 0.005; Chi2Cut : 100; MinClustersOnTrack : 4; MaxDistance : 0.02; HighPTCut : 10.0;
    @Flags : HighPTFit, VertexToTracker
    @Functions : CombineCollections, BuildNewTracks
    [VXDEncap]
    @Collections : VXDEndcapTrackerHits
    @Parameters : MaxCellAngle : 0.005; MaxCellAngleRZ : 0.005; Chi2Cut : 100; MinClustersOnTrack : 4; MaxDistance : 0.02; HighPTCut : 10.0;
    @Flags : HighPTFit, VertexToTracker
    @Functions : CombineCollections, ExtendTracks
    [LowerCellAngle1]
    @Collections : VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters : MaxCellAngle : 0.025; MaxCellAngleRZ : 0.025; Chi2Cut : 100; MinClustersOnTrack : 4; MaxDistance : 0.02; HighPTCut : 10.0;
    @Flags : HighPTFit, VertexToTracker, RadialSearch
    @Functions : CombineCollections, BuildNewTracks
    [LowerCellAngle2]
    @Collections :
    @Parameters : MaxCellAngle : 0.05; MaxCellAngleRZ : 0.05; Chi2Cut : 2000; MinClustersOnTrack : 4; MaxDistance : 0.02; HighPTCut : 10.0;
    @Flags : HighPTFit, VertexToTracker, RadialSearch
    @Functions : BuildNewTracks, SortTracks
    [Tracker]
    @Collections : ITrackerHits, OTrackerHits, ITrackerEndcapHits, OTrackerEndcapHits
    @Parameters : MaxCellAngle : 0.05; MaxCellAngleRZ : 0.05; Chi2Cut : 2000; MinClustersOnTrack : 4; MaxDistance : 0.02; HighPTCut : 1.0;
    @Flags : HighPTFit, VertexToTracker, RadialSearch
    @Functions : CombineCollections, ExtendTracks
    [Displaced]
    @Collections : VXDTrackerHits, VXDEndcapTrackerHits, ITrackerHits, OTrackerHits, ITrackerEndcapHits, OTrackerEndcapHits
    @Parameters : MaxCellAngle : 0.05; MaxCellAngleRZ : 0.05; Chi2Cut : 1000; MinClustersOnTrack : 5; MaxDistance : 0.015; HighPTCut : 10.0;
    @Flags : OnlyZSchi2cut, RadialSearch
    @Functions : CombineCollections, BuildNewTracks
)RAW"
  };
  // clang-format on
};

#endif  // ConformalTrackingV2_h
