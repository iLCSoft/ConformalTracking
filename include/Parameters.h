#ifndef STEPPARAMETERS_H
#define STEPPARAMETERS_H 1

#include <vector>

struct Parameters {
public:
  std::vector<int> _collections;  /// which collections to combine
  double           _maxCellAngle;
  double           _maxCellAngleRZ;
  double           _chi2cut;
  int              _minClustersOnTrack;
  double           _maxDistance;
  bool             _highPTfit;
  bool             _onlyZSchi2cut;
  bool             _radialSearch;
  bool             _vertexToTracker;
  int              _step;
  bool             _combine;
  bool             _build;
  bool             _extend;
  bool             _sortTracks;
  double           _tightenStep = 1;

  Parameters(std::vector<int> collections, double maxCellAngle, double maxCellAngleRZ, double chi2cut,
             int minClustersOnTrack, double maxDistance, bool highPTfit, bool onlyZSchi2cut, bool radialSearch,
             bool vertexToTracker, int step, bool combine, bool build, bool extend, bool sortTracks)
      : _collections(collections),
        _maxCellAngle(maxCellAngle),
        _maxCellAngleRZ(maxCellAngleRZ),
        _chi2cut(chi2cut),
        _minClustersOnTrack(minClustersOnTrack),
        _maxDistance(maxDistance),
        _highPTfit(highPTfit),
        _onlyZSchi2cut(onlyZSchi2cut),
        _radialSearch(radialSearch),
        _vertexToTracker(vertexToTracker),
        _step(step),
        _combine(combine),
        _build(build),
        _extend(extend),
        _sortTracks(sortTracks) {}

  Parameters(Parameters const& rhs) = default;

  void tighten() {
    double factor = (10.0 - _tightenStep) / (10.0 - (_tightenStep - 1.0));
    _maxCellAngle *= factor;
    _maxCellAngleRZ *= factor;
    _chi2cut *= factor;
    _tightenStep += 1.0;
  }
};

#endif  // STEPPARAMETERS_H
