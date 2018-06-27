#include "ParameterParser.h"

#include <marlin/Exceptions.h>
#include <streamlog/streamlog.h>

#include <sstream>

using PPVec     = std::vector<ParameterParser::ParsedParameters>;
using StringVev = std::vector<std::string>;

void ParameterParser::parseParameters(std::vector<Parameters>& parameters, StringVec const& rawSteps,
                                      StringVec const& allCollections) {
  std::stringstream configStream;
  for (auto const& line : rawSteps) {
    configStream << line << " ";
  }
  std::string configString = configStream.str();

  using It = std::string::const_iterator;
  ParameterGrammar<It> const myGrammar;

  PPVec parsedParameters;
  It    startString = std::begin(configString), endString = std::end(configString);
  bool  ok = qi::phrase_parse(startString, endString, myGrammar, qi::blank, parsedParameters);

  if (startString != endString) {
    std::stringstream errorMessage;
    errorMessage << "Remaining unparsed: ";
    while (startString != endString) {
      char c = *startString++;
      if (isprint(c)) {
        errorMessage << c;
      } else {
        errorMessage << "\\x" << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(c);
      }  //if
    }    //while
    streamlog_out(ERROR) << errorMessage.str();
    throw marlin::ParseException(errorMessage.str());
  }  //if startString!=endString

  if (not ok) {
    std::stringstream error;
    error << "Failed to parse step configuration: Check your config";
    throw marlin::ParseException(error.str());
  }  // not OK

  int count = -1;
  for (auto const& ps : parsedParameters) {
    parameters.emplace_back(ps, allCollections, ++count);
    streamlog_out(MESSAGE) << "Step " << count << ":" << ps << std::endl;
  }
  if (count == -1) {
    std::stringstream error;
    error << "Failed to parse step configuration: Check your config";
    throw marlin::ParseException(error.str());
  }

}  // parseParameters()
