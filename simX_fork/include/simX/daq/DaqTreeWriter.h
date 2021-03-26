//
// Created by munk on 23-09-15.
//

#ifndef SIMX_DAQTREEWRITER_H
#define SIMX_DAQTREEWRITER_H

#include "simX/Detection/Detector.h"
#include "simX/Detection/DetectionSystem.h"
#include "simX/parser/JsonOutputParser.h"
#include "DAQ.h"

#include <string>
#include <memory>
#include <Rtypes.h>
#include <ausa/TFileWrapper.h>

class TTree;
class TFile;

namespace simX {
    namespace detection {
        class Scaler;
    }

    namespace daq {
        class DaqTreeWriter {
        public:
            DaqTreeWriter(std::shared_ptr<AUSA::TFileWrapper> output,
                            const std::vector<std::shared_ptr<Detector>>& system,
                            const std::vector<std::shared_ptr<detection::Scaler>>& scalers,
                            const std::string& ausalibFile);

            DaqTreeWriter(std::shared_ptr<AUSA::TFileWrapper> output,
                          const detection::DetectionSystem& system,
                          const std::string& ausalibFile);

            ~DaqTreeWriter();

            void feed(DAQ::Output& output);
            void clear();


        private:
            std::vector<UInt_t> fMul, bMul;
            std::vector<std::vector<UInt_t>> fSeg, bSeg;
            std::vector<std::vector<UInt_t>> fEne, bEne;
            std::vector<std::vector<UInt_t>> fTdc, bTdc;

            std::vector<UInt_t> frontStrips;

            std::vector<UInt_t*> manualClear;

            std::shared_ptr<AUSA::TFileWrapper> wrapper;
            TFile* file;
            TTree* tree;

            const std::vector<std::shared_ptr<detection::Scaler>>& scalers;

            void map(UInt_t* p, const std::string& b, const std::string& mul);
            void map(const std::string& b, UInt_t* p);

            void mapMapping(parser::JsonOutputParser::Mapping& m, UInt_t* mul, UInt_t* seg, UInt_t* adc, UInt_t* tdc);
        };
    }
}
#endif //SIMX_DAQTREEWRITER_H
