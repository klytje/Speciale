//
// Created by munk on 20-11-15.
//

#ifndef SIMX_PARSER_IONPARSER_H
#define SIMX_PARSER_IONPARSER_H

#include "simX/parser/grammar/IonGrammar.h"

namespace simX {
    namespace parser {
        Ion parseIon(const std::string& input);
    }
}

#endif //SIMX_PARSER_IONPARSER_H
