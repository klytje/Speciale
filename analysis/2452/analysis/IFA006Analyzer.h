
#ifndef AUSAEVENTBUILDER_EXAMPLEANALYZER_H
#define AUSAEVENTBUILDER_EXAMPLEANALYZER_H

#include <TClonesArray.h>
#include "ausa/event/analyzer/Analyzer.h"
#include "ausa/util/OutputCollection.h"
#include "ausa/util/DynamicBranchVector.h"

class TTree;

namespace AUSA {
    namespace Event {

        /**
         * Example implementation of Analyzer - modify according to your needs!
         */
        class IFA006Analyzer : public Analyzer {

        public:
            IFA006Analyzer();
            virtual ~IFA006Analyzer();

            virtual void setup(std::shared_ptr<Setup> setup) override;
	        virtual void setScalers(const ScalerOutput&) override;
            virtual void analyze(const std::vector<PhysicsEvent> &events) override;

            virtual void terminate() override;

            virtual void saveToRootFile(TFileWrapper &file) override;

            virtual void reset() override;

        private:
            ScalerOutput scaler;
            size_t multCounter[6];
            AUSA::OutputCollection oc;
            std::vector<TH2*> h_loss;
            std::unique_ptr<TTree> tree;
            int dIndex[3], N, FI[3], BI[3], selected[3];
            double eCM[3], thetaCM[3], eLab[3], thetaLab[3], eDep[3], vZ[3], vZL[3], phiCM[3], phiLab[3], exBe8[3], prob[3], time[3];
            TClonesArray p, pCM;
            TLorentzVector pBe8;
            UInt_t mul;
            double exC12,  //exBe8
                    pTot, pX, pY, pZ, sumAng, devAng, deltaE, zangle, triggerw;
        };
    }
}

#endif //AUSAEVENTBUILDER_EXAMPLEANALYZER_H
