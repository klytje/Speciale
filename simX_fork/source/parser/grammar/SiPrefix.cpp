//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/SiPrefix.h"

simX::parser::SiPrefix::SiPrefix() {
    add
        ("Y", 24)
            ("Y", 24)
            ("Z", 21)
            ("E", 18)
            ("P", 15)
            ("T", 12)
            ("G", 9)
            ("M", 6)
            ("k", 3)
            ("h", 2)
            ("da", 1)
            ("d", -1)
            ("c", -2)
            ("m", -3)
            ("u", -6)
            ("n", -9)
            ("p", -12)
            ("f", -15)
            ("a", -18)
            ("z", -21)
            ("y", -24)
        ;

    name("SI prefix");
}
