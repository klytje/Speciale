//
// Created by jesper on 11/5/15.
//

#include "simX/Detection/Detector.h"

std::string simX::Detector::description() {
    return getName();
}


namespace simX{
    std::ostream& operator<<(std::ostream& stream, Detector& det) {
        stream << det.description();
        return stream;
    }
}
