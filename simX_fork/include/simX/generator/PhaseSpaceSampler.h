#ifndef SIMX_PHASE_SPACE_SAMPLER_H
#define	SIMX_PHASE_SPACE_SAMPLER_H

#include "FinalStateGenerator.h"
#include "simX/parser/ReactionParser.h"
#include <memory>

namespace simX {
    class Particle;

    namespace samplers {
        class SampleSet;
    }

    class PhaseSpaceSampler : public FinalStateGenerator {
    public:
        using Input = std::unique_ptr<samplers::SampleSet>;
        PhaseSpaceSampler(Input input_, std::vector<Particle*> daughters, std::string unit);

        virtual std::vector<TLorentzVector> &getFourMomenta() override;

        static std::unique_ptr<FinalStateGenerator> build(const NBodyDecay& d, parser::ReactionParser::Options& opt);

    private:
        Input input;
        const std::vector<Particle*> daughters;
        std::vector<TLorentzVector> result;
        double multiplier;
    };
}

#endif // SIMX_PHASE_SPACE_SAMPLER_H