//
// Created by munk on 29-06-15.
//

#ifndef SIMX_KEYVALUEGRAMMAR_H
#define SIMX_KEYVALUEGRAMMAR_H

#include <boost/spirit/include/qi.hpp>
#include <map>
#include <string>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;
        namespace ascii = boost::spirit::ascii;

        typedef std::string::const_iterator Iterator;
        typedef std::map<std::string, std::string> KeyValueMap;

        struct KeyValueGrammar : qi::grammar<Iterator, KeyValueMap(), ascii::space_type>
        {
            KeyValueGrammar();

            qi::rule<Iterator, KeyValueMap(), ascii::space_type> start;
            qi::rule<Iterator, std::pair<std::string, std::string>(), ascii::space_type> pair;
            qi::rule<Iterator, std::string(), ascii::space_type> key, value;
        };
    }
}
#endif //SIMX_KEYVALUEGRAMMAR_H
