//
// Created by jesper on 4/4/16.
//

#ifndef SIMX_SINGLEBODYKINEMATICS_H
#define SIMX_SINGLEBODYKINEMATICS_H

#include "simX/generator/FinalStateGenerator.h"
#include "simX/Particle.h"

namespace simX {

    /**
    * Dummy implementation of a single body reaction, i.e. it simply forwards the momenta of the incomming particle
    * two-body kinematics.
    * @param particle Particle to be forwarded.
    */
    class SingleBodyKinematics : public FinalStateGenerator {

    public:
        SingleBodyKinematics( const Particle &particle);
        virtual ~SingleBodyKinematics() = default;

        virtual std::vector<TLorentzVector>& getFourMomenta() override;

    private:
        const Particle& particle;
        std::vector<TLorentzVector> momenta;
    };
}

#endif //SIMX_SINGLEBODYKINEMATICS_H
