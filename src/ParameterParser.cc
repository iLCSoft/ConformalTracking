#include "ParameterParser.h"

#include <marlin/Exceptions.h>
#include <streamlog/streamlog.h>

#include <sstream>
#include <stdexcept>

// we put the grammar into the implementation to avoid excessive compilation time

using PPVec     = std::vector<ParameterParser::ParsedParameters>;
using StringVev = std::vector<std::string>;

namespace qi    = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

//Create a grammar that returns a vector of ParsedParameters
template <typename It, typename Skipper = qi::blank_type>
struct ParameterGrammar : qi::grammar<It, ParameterParser::PPVec(), Skipper> {
  ParameterGrammar() : ParameterGrammar::base_type(start) {
    // name is a sequence of characters/numbers/underscore and hyphen
    name = +qi::char_("-a-zA-Z_0-9");

    //separation can be done by arbirary number of newline, semicolon, comma, or colon
    spacer = *(qi::eol | ';' | ',' | ':');

    //the name of the step is a name in brackets
    stepName = qi::lit("[") >> name >> qi::lit("]");

    //collections are a list of names separated by spacers, can be empty
    collections = *(name % spacer);
    //functions and flags are names separated by spacers
    functions = name % spacer;
    flags     = *(name % spacer);

    //a key cannot start with a number
    key = qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");

    // a parameters is a key followed by equal sign and a double
    parPair = key >> -('=' >> qi::double_);
    //parameters is a parsing pairs into a map, separated by spacers
    parameters = parPair % spacer;

    // clang-format off
    // The order is the same as the struct we pass into its, a number of ParsedParameters
    // described by a ordered! sequence of Collections, parameters, flags and functions
    start = *(*qi::eol >> stepName
              >> spacer >> "@Collections" >> qi::lit('=') >> collections
              >> spacer >> "@Parameters" >> qi::lit('=') >> parameters
              >> spacer >> "@Flags" >> qi::lit('=') >> flags
              >> spacer >> "@Functions" >> qi::lit('=') >> functions
              >> spacer);
    // clang-format on

    // Uncomment to see spirit at work

    //BOOST_SPIRIT_DEBUG_NODES((start)(collections)(flags)(parameters)(functions)(parPair)(key)(name)(stepName));
  }

private:
  qi::rule<It, PPVec(), Skipper>                 start{};
  qi::rule<It, Parameters::StringVec(), Skipper> collections{};
  qi::rule<It, Parameters::ParMap(), Skipper>    parameters{};
  qi::rule<It, std::string(), Skipper>           name{};
  qi::rule<It, std::string(), Skipper>           stepName{};
  qi::rule<It, std::string(), Skipper>           key{};
  qi::rule<It, std::pair<std::string, double>(), Skipper> parPair{};
  qi::rule<It, Parameters::StringVec(), Skipper> flags{};
  qi::rule<It, Parameters::StringVec(), Skipper> functions{};
  qi::rule<It, Skipper> spacer{};
};

void ParameterParser::parseParameters(std::vector<Parameters>& parameters, StringVec const& rawSteps,
                                      StringVec const& allCollections) {
  //concatenate the StringVec into a single string
  std::stringstream configStream;
  for (auto const& line : rawSteps) {
    configStream << line << " ";
  }
  std::string configString = configStream.str();

  //Parse the string into the ParsedParameterVector
  PPVec parsedParameters;
  using It = std::string::const_iterator;
  ParameterGrammar<It> const myGrammar;
  It                         startString = std::begin(configString), endString = std::end(configString);
  bool                       ok = qi::phrase_parse(startString, endString, myGrammar, qi::blank, parsedParameters);

  //if there was an error during parsing
  if (not ok) {
    std::stringstream error;
    error << "Failed to parse step configuration: Check your config";
    throw marlin::ParseException(error.str());
  }  // not OK

  //Make sure everything was parsed, if not throw
  if (startString != endString) {
    std::stringstream errorMessage;
    errorMessage << "Error parsing parameters: remaining unparsed:\n";
    while (startString != endString) {
      errorMessage << *startString++;
    }
    streamlog_out(ERROR) << errorMessage.str();
    throw marlin::ParseException(errorMessage.str());
  }  //if startString!=endString

  //Transfer the ParsedParameters into the Parameters
  int count = -1;
  for (auto const& ps : parsedParameters) {
    try {
      parameters.emplace_back(ps, allCollections, ++count);
      streamlog_out(DEBUG8) << "Step " << count << ":" << ps << std::endl;
    } catch (std::out_of_range& e) {
      std::stringstream error;
      error << "Failed to create step configuration: A parameter is missing for " << ps._name;
      throw marlin::ParseException(error.str());
    }
  }

  //if there were no parameters in the config
  if (count == -1) {
    std::stringstream error;
    error << "Failed to parse step configuration: Check your config";
    throw marlin::ParseException(error.str());
  }

}  // parseParameters()
