#ifndef PARAMETERPARSER_H
#define PARAMETERPARSER_H 1

#include "Parameters.h"

//#define BOOST_SPIRIT_DEBUG

// Prevent inclusing of boost/phoenix/stl/tuple.hpp
// which is causing issues like
// multiple definition of `boost::phoenix::placeholders::uarg3'; CMakeFiles/ConformalTracking.dir/ConformalTracking/src/ConformalTrackingV2.cc.o:(.bss+0x636): first defined here

#define BOOST_PHOENIX_STL_TUPLE_H_

#include <boost/fusion/include/std_pair.hpp>
#include <boost/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace ParameterParser {

  struct ParsedParameters {
  public:
    std::string           _name{};
    Parameters::StringVec _collections{};  /// which collections to combine
    Parameters::ParMap    _parameters{};
    Parameters::StringVec _flags{};
    Parameters::StringVec _functions{};

    friend std::ostream& operator<<(std::ostream& out, ParsedParameters const& p) {
      out << "Parameters " << p._name << ":\n";
      out << "\tCollections:";
      for (auto const& col : p._collections) {
        out << "\n\t\t" << col;
      }
      out << "\n\tParameters:";
      for (auto const& par : p._parameters) {
        out << "\n\t\t" << par.first << ":" << par.second;
      }
      out << "\n\tFlags: ";
      for (auto const& flag : p._flags) {
        out << "\n\t\t" << flag;
      }
      out << "\n\tFunctions:";
      for (auto const& func : p._functions) {
        out << "\n\t\t" << func;
      }
      return out;
    }
  };
  using PPVec = std::vector<ParsedParameters>;
}  // namespace ParameterParser

BOOST_FUSION_ADAPT_STRUCT(ParameterParser::ParsedParameters,
                          (std::string, _name)(Parameters::StringVec, _collections)(Parameters::ParMap, _parameters)(
                              Parameters::StringVec, _flags)(Parameters::StringVec, _functions))

namespace ParameterParser {

  void parseParameters(std::vector<Parameters>& parameters, std::vector<std::string> const& rawSteps,
                       std::vector<std::string> const& allCollections);

}  // namespace ParameterParser
#endif  // PARAMETERPARSER_H
