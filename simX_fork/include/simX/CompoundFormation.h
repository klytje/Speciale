//
// Created by munk on 28-06-15.
//

#ifndef SIMX_FUSIONREACTION_H
#define SIMX_FUSIONREACTION_H

#include "simX/NuclearProcess.h"
#include "simX/weight/WeightCalculator.h"
#include "simX/generator/FinalStateGenerator.h"

#include <memory>

namespace simX {
    class CompoundFormation : public NuclearProcess {
    public:
        CompoundFormation(Particle& beam, Particle target);

        virtual void runProcess() override;
        virtual std::vector<Particle *>& getDaughters() override;
        virtual const std::vector<Particle *>& getDaughters() const override;
        virtual double getWeight() const override;
        virtual void setWeightCalculator( std::unique_ptr<WeightCalculator> ) override;
        virtual int getNumberOfNucleons() const override;
        virtual int getNumberOfProtonsInInitialState() const override;        

    private:
        std::vector<Particle *> daughters;
        Particle target, product;
        const Particle& beam;
        std::unique_ptr<WeightCalculator> weightCalculator;
        std::unique_ptr<FinalStateGenerator> fsGenerator;
    };
}
#endif //SIMX_FUSIONREACTION_H
