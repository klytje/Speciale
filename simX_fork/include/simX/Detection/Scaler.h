//
// Created by munk on 17-02-16.
//

#ifndef SIMX_SCALER_H
#define SIMX_SCALER_H

#include <string>
#include <Rtypes.h>

namespace simX {
    namespace detection {
        class Scaler {
        public:
            Scaler(const std::string& name, UInt_t v = 0) : name(name), value(v) {

            }

            const std::string& getName() const {
                return name;
            }

            UInt_t getValue() const {
                return value;
            }

            UInt_t* getPointer() {
                return &value;
            }

            const UInt_t* getPointer() const {
                return &value;
            }

        private:
            std::string name;
            UInt_t value;
        };
    }
}
#endif //SIMX_SCALER_H
