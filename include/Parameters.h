#ifndef STEPPARAMETERS_H
#define STEPPARAMETERS_H 1

#include <map>
#include <string>
#include <vector>

namespace ParameterParser {
  struct ParsedParameters;
}

struct Parameters {
  using ParMap    = std::map<std::string, double>;
  using StringVec = std::vector<std::string>;

public:
  std::vector<int> _collections;  /// which collections to combine
  double           _maxCellAngle;
  double           _maxCellAngleRZ;
  double           _chi2cut;
  int              _minClustersOnTrack;
  double           _maxDistance;
  double           _highPTcut;
  bool             _highPTfit;
  bool             _onlyZSchi2cut;
  bool             _radialSearch;
  bool             _vertexToTracker;
  bool             _kalmanFitForward;
  int              _step;
  bool             _combine;
  bool             _build;
  bool             _extend;
  bool             _sortTracks;
  double           _tightenStep = 1;

  const StringVec _existingFunctions = {
      "CombineCollections", "ExtendTracks", "BuildNewTracks", "SortTracks",
  };
  const StringVec _existingFlags = {
      "HighPTFit", "OnlyZSchi2cut", "RadialSearch", "VertexToTracker", "KalmanFitForward", "KalmanFitBackward",
  };
  const StringVec _existingParameters = {
      "MaxCellAngle", "MaxCellAngleRZ", "Chi2Cut", "MinClustersOnTrack", "MaxDistance", "HighPTCut",
  };

  Parameters(std::vector<int> const& collections, double maxCellAngle, double maxCellAngleRZ, double chi2cut,
             int minClustersOnTrack, double maxDistance, double highPTcut, bool highPTfit, bool onlyZSchi2cut,
             bool radialSearch, bool vertexToTracker, bool kalmanFitForward, int step, bool combine, bool build, bool extend,
             bool sortTracks)
      : _collections(collections),
        _maxCellAngle(maxCellAngle),
        _maxCellAngleRZ(maxCellAngleRZ),
        _chi2cut(chi2cut),
        _minClustersOnTrack(minClustersOnTrack),
        _maxDistance(maxDistance),
        _highPTcut(highPTcut),
        _highPTfit(highPTfit),
        _onlyZSchi2cut(onlyZSchi2cut),
        _radialSearch(radialSearch),
        _vertexToTracker(vertexToTracker),
        _kalmanFitForward(kalmanFitForward),
        _step(step),
        _combine(combine),
        _build(build),
        _extend(extend),
        _sortTracks(sortTracks) {}

  explicit Parameters(ParameterParser::ParsedParameters const& ps, std::vector<std::string> const& allCollections,
                      int step = 0);

  Parameters(Parameters const& rhs) = default;

  void tighten() {
    double factor = (10.0 - _tightenStep) / (10.0 - (_tightenStep - 1.0));
    _maxCellAngle *= factor;
    _maxCellAngleRZ *= factor;
    _chi2cut *= factor;
    _tightenStep += 1.0;
  }

private:
  void check(StringVec const& values, StringVec const& options, std::string const& type);
  void check(ParMap const& values, StringVec const& options, std::string const& type);
};

#endif  // STEPPARAMETERS_H
