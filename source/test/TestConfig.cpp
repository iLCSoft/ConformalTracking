#define BOOST_SPIRIT_DEBUG 1
#include "ParameterParser.h"

using PPVec  = std::vector<ParameterParser::ParsedParameters>;
namespace qi = boost::spirit::qi;

//Create a grammar that returns a vector of ParsedParameters
template <typename It, typename Skipper = qi::blank_type>
struct ParameterGrammar : qi::grammar<It, ParameterParser::PPVec(), Skipper> {
  ParameterGrammar() : ParameterGrammar::base_type(start) {
    // name is a sequence of characters/numbers/underscore and hyphen
    name = +qi::char_("-a-zA-Z_0-9");

    //separation can be done by arbirary number of newline, semicolon, comma, or colon
    spacer = *(qi::eol | ';' | ',' | ':');

    //assignment can be done with = or colon
    equals = (qi::lit('=') | qi::lit(':'));

    //the name of the step is a name in brackets
    stepName = qi::lit("[") >> name >> qi::lit("]");

    //collections are a list of names separated by spacers, can be empty
    collections = *(name % spacer);
    //functions and flags are names separated by spacers
    functions = name % spacer;
    flags     = *(name % spacer);

    //a key cannot start with a number
    key = qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");

    // a parameters is a key followed by equals (=|:) and a double
    parPair = key >> -(equals >> qi::double_);
    //parameters is a parsing pairs into a map, separated by spacers
    parameters = parPair % spacer;

    // clang-format off
    // The order is the same as the struct we pass into its, a number of ParsedParameters
    // described by a ordered! sequence of Collections, parameters, flags and functions
    start = *(*qi::eol >> stepName
              >> spacer >> "@Collections" >> equals >> collections
              >> spacer >> "@Parameters" >> equals >> parameters
              >> spacer >> "@Flags" >> equals >> flags
              >> spacer >> "@Functions" >> equals >> functions
              >> spacer);
    // clang-format on

    // Uncomment to see spirit at work

    //BOOST_SPIRIT_DEBUG_NODES((start)(collections)(flags)(parameters)(functions)(parPair)(key)(name)(stepName));
  }

private:
  qi::rule<It, PPVec(), Skipper>                          start{};
  qi::rule<It, Parameters::StringVec(), Skipper>          collections{}, flags{}, functions{};
  qi::rule<It, Parameters::ParMap(), Skipper>             parameters{};
  qi::rule<It, std::string(), Skipper>                    name{}, stepName{}, key{};
  qi::rule<It, std::pair<std::string, double>(), Skipper> parPair{};
  qi::rule<It, Skipper>                                   spacer{}, equals{};
};

int main() {
  using It = std::string::const_iterator;
  ParameterGrammar<It> const myGrammar;

  std::string configString = R"RAW(
    [VXDBarrel]
    @Collections : VXDTrackerHits
    @Parameters : MaxCellAngle : 0.035; MaxCellAngleRZ : 0.035;  Chi2Cut : 300;  MinClustersOnTrack : 5;  MaxDistance : 0.015
    @Flags : HighPTFit, VertexToTracker, KalmanFitForward
    @Functions : CombineCollections, BuildNewTracks
    [VXDEncap]
    @Collections : VXDEndcapTrackerHits
    @Parameters : MaxCellAngle : 0.035 MaxCellAngleRZ : 0.035:  Chi2Cut : 300:  MinClustersOnTrack : 5   MaxDistance : 0.015
    @Flags : HighPTFit, VertexToTracker, KalmanFitForward
    @Functions : CombineCollections, ExtendTracks
    [TCVC]
    @Collections : VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters : MaxCellAngle : 0.035; MaxCellAngleRZ : 0.035;  Chi2Cut : 300;  MinClustersOnTrack : 5;  MaxDistance : 0.015
    @Flags : HighPTFit, VertexToTracker, KalmanFitForward
    @Functions : CombineCollections, BuildNewTracks
    [Foo]
    @Collections :
    @Parameters : MaxCellAngle : 0.035; MaxCellAngleRZ : 0.035;  Chi2Cut : 300;  MinClustersOnTrack : 5;  MaxDistance : 0.015
    @Flags : HighPTFit, VertexToTracker, KalmanFitBackward
    @Functions : BuildNewTracks
)RAW";

  ParameterParser::PPVec parsed;
  It                     startString = begin(configString), endString = end(configString);
  bool                   ok = qi::phrase_parse(startString, endString, myGrammar, qi::blank, parsed);

  std::cout << "--- File ---"
            << "\n";
  std::cout << "parsed   = " << std::boolalpha << ok << "\n";
  if (ok) {
    std::cout << "Number " << parsed.size() << std::endl;
    int count = 0;
    for (auto const& step : parsed) {
      std::cout << "Step" << count++ << ":" << step << std::endl;
    }
  } else {
    return 1;
  }

  if (startString != endString) {
    std::cout << "Remaining unparsed: ";
    while (startString != endString) {
      char c = *startString++;
      if (isprint(c)) {
        std::cout << c;
      } else {
        std::cout << "\\x" << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(c);
      }  //if
    }  //while
    std::cout << std::endl;
    return 1;
  }  //if
  std::cout << std::endl;
  return 0;
}  //main
