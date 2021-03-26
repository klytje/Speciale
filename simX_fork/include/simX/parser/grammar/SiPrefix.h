//
// Created by munk on 25-06-15.
//

#ifndef SIMX_SIPREFIX_H
#define SIMX_SIPREFIX_H

#include <boost/spirit/include/qi.hpp>

namespace simX {
    namespace parser {
        namespace qi = boost::spirit::qi;

        struct SiPrefix : qi::symbols<char, int> {
            SiPrefix();
        };
    }
}
#endif //SIMX_SIPREFIX_H
