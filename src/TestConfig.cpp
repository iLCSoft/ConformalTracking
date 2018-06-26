#define BOOST_SPIRIT_DEBUG 1
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <iomanip>
#include <boost/fusion/include/std_pair.hpp>

#include <algorithm>
#include <vector>
#include <map>
#include <stdexcept>


namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

using ParMap = std::map<std::string, double>;
using FlagMap = std::map<std::string, bool>;

struct Parameters_t {
public:
  std::vector<std::string> _collections{};  /// which collections to combine
  ParMap _parameters{};
  FlagMap _flags{};
  std::vector<std::string> _functions{};

  static void allowedFunctions(std::string const& func){
    static std::vector<std::string> af = {"CombineCollections",
                                          "ExtendTracks",
                                          "BuildNewTracks",
    };
    if(std::find(af.begin(), af.end(), func) == af.end()) {
      throw std::runtime_error("Function not allowed");
    }
  }

  friend  std::ostream& operator<<(std::ostream& out, Parameters_t const&p){
    std::cout << "PARAMETERS:\n";
    out << "\tCollections:";
    for (auto const& col: p._collections) {
      out << "\n\t\t" << col;
    }
    out << "\n\tParameters:";
    for (auto const& par: p._parameters) {
      out << "\n\t\t" << par.first << ":" << par.second;
    }
    out << "\n\tFlags:";
    for (auto const& flag: p._flags) {
      out << "\n\t\t" << flag.first << ":" << flag.second;
    }
    out << "\n\tFunctions:";
    for (auto const& func: p._functions) {
      out << "\n\t\t" << func;
    }
    return out;
  }

};


BOOST_FUSION_ADAPT_STRUCT(Parameters_t,
                          (std::vector<std::string>, _collections )
                          (ParMap, _parameters )
                          (FlagMap, _flags)
                          (std::vector<std::string>, _functions )
                          )

template <typename It, typename Skipper = qi::blank_type>
struct grammar : qi::grammar<It, std::vector<Parameters_t>(), Skipper> {

  grammar() : grammar::base_type(start) {

    name = +qi::char_("-a-zA-Z_0-9");
    collections = name % ',';
    functions = name % ',';

    keyRule = qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");

    parameters = parPair % (qi::lit(';') | qi::lit(','));
    parPair = keyRule >> -('=' >> qi::double_);

    flags = flagPair % (qi::lit(';') | qi::lit(','));
    flagPair = keyRule >> -('=' >> qi::bool_);

    start = *(qi::no_case["[STEP]"]
              >> qi::eol >> "Collections" >> qi::lit('=') >> collections
              >> qi::eol >> "Parameters" >> qi::lit('=') >> parameters
              >> qi::eol >> "Flags" >> qi::lit('=') >> flags
              >> qi::eol >> "Functions" >> qi::lit('=') >> functions
              >> *(qi::char_(';') | qi::eol))
    ;

    BOOST_SPIRIT_DEBUG_NODES((start)(collections)(parameters)(functions)(parPair)(keyRule)(name));

  }

private:
  qi::rule<It, std::vector<Parameters_t>(), Skipper> start{};
  qi::rule<It, std::vector<std::string>(), Skipper> collections{};
  qi::rule<It, ParMap(), Skipper> parameters{};
  qi::rule<It, std::string(), Skipper> name{};
  qi::rule<It, std::string(), Skipper> keyRule{};
  qi::rule<It, std::pair<std::string, double>(), Skipper> parPair{};
  qi::rule<It, FlagMap(), Skipper> flags{};
  qi::rule<It, std::pair<std::string, bool>(), Skipper> flagPair{};
  qi::rule<It, std::vector<std::string>(), Skipper> functions{};


};

int main() {
  using It = std::string::const_iterator;
  grammar<It> const myGrammar;

  std::string configString = R"RAW([Step]
Collections = Foo,Bar
Parameters= MaxAngle = 2000;MaxAngleRZ = -50;Chi2Cut = 80;MinCluOnTrack=5;MaxDist=123
Flags=HighPTFit=false;OnlyZSChi2=false;RadialSearch=true;VertexToTracker=false
Functions = CombineCollections,BuildNewTracks
[Step]
Collections = Baz
Parameters  = MaxAngle = 2433;MaxAngleRZ = -2350;Chi2Cut = 1280;MinCluOnTrack=12;MaxDist=123
Flags       = HighPTFit=false;OnlyZSChi2=false;RadialSearch=true;VertexToTracker=false
Functions   = CombineCollections,ExtendTracks
)RAW";


  std::vector<Parameters_t> parsed;
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
