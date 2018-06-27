#ifndef PARAMETERPARSER_H
#define PARAMETERPARSER_H 1

#include "Parameters.h"

//#define BOOST_SPIRIT_DEBUG
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/phoenix.hpp>
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
}  //namespace

BOOST_FUSION_ADAPT_STRUCT(ParameterParser::ParsedParameters,
                          (std::string, _name)(Parameters::StringVec, _collections)(Parameters::ParMap, _parameters)(
                              Parameters::StringVec, _flags)(Parameters::StringVec, _functions))

namespace ParameterParser {

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

  void parseParameters(std::vector<Parameters>& parameters, std::vector<std::string> const& rawSteps,
                       std::vector<std::string> const& allCollections);

}  //namespace
#endif  // PARAMETERPARSER_H
