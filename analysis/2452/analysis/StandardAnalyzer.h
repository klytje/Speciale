
#ifndef AUSAEVENTBUILDER_EXAMPLEANALYZER_H
#define AUSAEVENTBUILDER_EXAMPLEANALYZER_H

#include "ausa/event/analyzer/Analyzer.h"

class TTree;

namespace AUSA {
    namespace Event {

        /**
         * Example implementation of Analyzer - modify according to your needs!
         */
        class StandardAnalyzer : public Analyzer {

        public:
            StandardAnalyzer();
            virtual ~StandardAnalyzer();

            virtual void setup(std::shared_ptr<Setup> setup) override;
	        virtual void setScalers(const ScalerOutput&) override;
            virtual void analyze(const std::vector<PhysicsEvent> &events) override;

            virtual void terminate() override;

            virtual void saveToRootFile(TFileWrapper &file) override;

            virtual void reset() override;

        private:
            size_t multCounter[6];
            int dIndex[3], FI[3], BI[3], selected[3];
            double eCM[3], thetaCM[3], eLab[3], thetaLab[3], eDep[3], vZ[3], vZL[3], phiCM[3], phiLab[3], time[3];
            TLorentzVector pBe8;
            double pTot, pX, pY, pZ, deltaE, triggerw;

            std::unique_ptr<TTree> tree;
            ScalerOutput scaler;
            int d[3], N, mul;
            double prob[3],
                E_cm[3], E_lab[3], E_dep[3], dE,
                vzl[3], vz[3], 
                theta_cm[3], theta_lab[3], 
                phi_cm[3], phi_lab[3],
                exBe8[3], exC12,
                px, py, pz, ptot,
                zangle, sumAng, devAng;
        };
    }
}

#endif //AUSAEVENTBUILDER_EXAMPLEANALYZER_H
