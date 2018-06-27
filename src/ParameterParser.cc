#include "ParameterParser.h"

#include <marlin/Exceptions.h>
#include <streamlog/streamlog.h>

#include <sstream>

using PPVec     = std::vector<ParameterParser::ParsedParameters>;
using StringVev = std::vector<std::string>;

namespace qi    = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template <typename It, typename Skipper = qi::blank_type>
struct ParameterGrammar : qi::grammar<It, ParameterParser::PPVec(), Skipper> {
  ParameterGrammar() : ParameterGrammar::base_type(start) {
    name = +qi::char_("-a-zA-Z_0-9");

    stepName = qi::lit("[") >> name >> qi::lit("]");

    collections = name % spacer;
    functions   = name % spacer;
    flags       = name % spacer;

    key = qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");

    parameters = parPair % (qi::lit(';') | qi::lit(','));
    parPair    = key >> -('=' >> qi::double_);

    spacer = *(qi::eol | ';' | ',' | ':');

    // clang-format off
    start = *(*qi::eol >> stepName
              >> spacer >> "@Collections" >> qi::lit('=') >> collections
              >> spacer >> "@Parameters" >> qi::lit('=') >> parameters
              >> spacer >> "@Flags" >> qi::lit('=') >> flags
              >> spacer >> "@Functions" >> qi::lit('=') >> functions
              >> spacer);
    // clang-format on

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
