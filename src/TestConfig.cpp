#define BOOST_SPIRIT_DEBUG 1
#include "ParameterParser.h"

int main() {
  using It = std::string::const_iterator;
  ParameterParser::ParameterGrammar<It> const myGrammar;

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

  namespace qi = ParameterParser::qi;
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
