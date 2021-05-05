
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
            std::unique_ptr<TTree> tree;
            ScalerOutput scaler;
            int mi[3], mul, N;
            double prob[3], 
                E_cm[3], E_lab[3], E_dep[3], dE,
                vzl[3], vz[3], 
                theta_cm[3], theta_lab[3], 
                phi_cm[3], phi_lab[3],
                exBe8[3], exC12,
                px[3], py[3], pz[3], ptot,
                zangle, sumAng, devAng;
        };
    }
}

#endif //AUSAEVENTBUILDER_EXAMPLEANALYZER_H
