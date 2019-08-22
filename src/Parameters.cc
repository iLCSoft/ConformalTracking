#include "Parameters.h"
#include "ParameterParser.h"

#include <marlin/Exceptions.h>

bool findEntry(std::vector<std::string> const& functions, std::string const& entry) {
  return std::find(functions.begin(), functions.end(), entry) != functions.end();
}

Parameters::Parameters(ParameterParser::ParsedParameters const& ps, std::vector<std::string> const& allCollections, int step)
    : _collections({}),
      _maxCellAngle(ps._parameters.at("MaxCellAngle")),
      _maxCellAngleRZ(ps._parameters.at("MaxCellAngleRZ")),
      _chi2cut(ps._parameters.at("Chi2Cut")),
      _minClustersOnTrack(ps._parameters.at("MinClustersOnTrack")),
      _maxDistance(ps._parameters.at("MaxDistance")),
      _maxSlopeZ(ps._parameters.at("SlopeZRange")),
      _highPTcut(ps._parameters.at("HighPTCut")),
      _highPTfit(findEntry(ps._flags, "HighPTFit")),
      _onlyZSchi2cut(findEntry(ps._flags, "OnlyZSchi2cut")),
      _radialSearch(findEntry(ps._flags, "RadialSearch")),
      _vertexToTracker(findEntry(ps._flags, "VertexToTracker")),
      _kalmanFitForward(findEntry(ps._flags, "KalmanFitForward") or not findEntry(ps._flags, "KalmanFitBackward")),
      _step(step),
      _combine(findEntry(ps._functions, "CombineCollections")),
      _build(findEntry(ps._functions, "BuildNewTracks")),
      _extend(findEntry(ps._functions, "ExtendTracks")),
      _sortTracks(findEntry(ps._functions, "SortTracks")) {
  for (auto const& colName : ps._collections) {
    auto it = std::find(allCollections.begin(), allCollections.end(), colName);
    if (it == allCollections.end()) {
      std::stringstream error;
      error << ": Collection name \"" << colName << "\" not found in Collection parameter! Check your config";
      throw marlin::ParseException(error.str());
    }
    _collections.push_back(it - allCollections.begin());

    check(ps._parameters, _existingParameters, "Parameter");
    check(ps._functions, _existingFunctions, "Function");
    check(ps._flags, _existingFlags, "Flag");
  }
}

void Parameters::check(Parameters::StringVec const& values, Parameters::StringVec const& options, std::string const& type) {
  for (auto const& value : values) {
    if (not findEntry(options, value)) {
      std::stringstream error;
      error << "Unknown " << type << ": " << value << ". Check your configuration";
      throw marlin::ParseException(error.str());
    }
  }
}

void Parameters::check(Parameters::ParMap const& values, Parameters::StringVec const& options, std::string const& type) {
  for (auto const& value : values) {
    if (not findEntry(options, value.first)) {
      std::stringstream error;
      error << "Unknown " << type << ": " << value.first << ". Check your configuration";
      throw marlin::ParseException(error.str());
    }
  }
}
