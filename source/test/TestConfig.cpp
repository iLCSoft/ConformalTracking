#define BOOST_SPIRIT_DEBUG 1
#include "ParameterParser.h"

namespace qi = boost::spirit::qi;

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
  qi::rule<It, ParameterParser::PPVec(), Skipper>                 start{};
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


int main() {
  using It = std::string::const_iterator;
  ParameterGrammar<It> const myGrammar;

  std::string configString = R"RAW(
    [VXDBarrel]
    @Collections = VXDTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, BuildNewTracks
    [VXDEncap]
    @Collections = VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, ExtendTracks
    [TCVC]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 5;  MaxDistance = 0.015
    @Flags = HighPTFit, VertexToTracker;
    @Functions = CombineCollections, BuildNewTracks
)RAW";

  ParameterParser::PPVec parsed;
  It startString = begin(configString), endString = end(configString);
  bool ok = qi::phrase_parse(startString, endString, myGrammar, qi::blank, parsed);

  std::cout << "--- File ---" << "\n";
  std::cout << "parsed   = " << std::boolalpha << ok << "\n";
  if (ok) {
    std::cout << "Number " << parsed.size()  << std::endl;
    int count = 0;
    for (auto const& step : parsed ) {
      std::cout << "Step" << count++ << ":" << step  << std::endl;
    }
  }

  if (startString!=endString) {
    std::cout << "Remaining unparsed: ";
    while (startString!=endString) {
      char c = *startString++;
      if (isprint(c)) {
        std::cout << c;
      } else {
        std::cout << "\\x" << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(c);
      } //if
    }//while
  }//if
    std::cout << std::endl;
    //}
}//main
