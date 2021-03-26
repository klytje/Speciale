
#ifndef NBODYDECAY_H
#define	NBODYDECAY_H

#include "simX/NuclearProcess.h"
#include "simX/Particle.h"
#include "simX/angular/FixedAngularCorrelation.h"

#include <memory>
#include "ausa/util/memory"



namespace simX {
    class FinalStateGenerator;
    class AmpltudeCalculator;

    /**
    * Derived class for simulating N-body breakups.
    * @param parent Reference to parent nucleus.
    * @param daughters Decay products.
    * @param ac Angular correlation function.
    * @param wCalc Weight calculator.
    */
    class NBodyDecay : public NuclearProcess {
        enum class Z_AXIS {
            COLINEAR,
            ORTHOGONAL
        };
    

        public:
            NBodyDecay( Particle& parent, std::vector<Particle> theDaughters,
                        std::unique_ptr<FinalStateGenerator> generator = nullptr,
                        std::unique_ptr<angular::AngularCorrelation> ac=std::make_unique<angular::FixedAngularCorrelation>(0,0),
                        std::unique_ptr<WeightCalculator> ampCalc = nullptr,
                        int orbitalL = 0, Z_AXIS axis = Z_AXIS::COLINEAR);

            static NBodyDecay withDefaultGenerator(Particle& parent, std::vector<Particle> theDaughters,
                                                   std::unique_ptr<angular::AngularCorrelation> ac=std::make_unique<angular::FixedAngularCorrelation>(0,0),
                                                   std::unique_ptr<WeightCalculator> ampCalc = nullptr);

            NBodyDecay(NBodyDecay&&);

            virtual ~NBodyDecay();
    
            virtual void runProcess() override;
            virtual std::vector<Particle*>& getDaughters() override;
            virtual const std::vector<Particle*>& getDaughters() const override;
            virtual double getWeight() const override;
            virtual void setWeightCalculator(std::unique_ptr<WeightCalculator>) override;
            virtual void setFinalStateGenerator( std::unique_ptr<FinalStateGenerator> generator);
            virtual int getNumberOfNucleons() const override;
            virtual int getNumberOfProtonsInInitialState() const override; 
        
            const WeightCalculator* const getWeightCalculator() const {
                return weightCalculator.get();
            }

            /**
             * Return the multiplicity of the decay (number of daughters)
             */
            inline int getMultiplicity() const { return MULT; };

            virtual int getL() const {return orbitalL;}
            
        protected:

        private:
            int MULT;
            double mDiff;
            const Particle& parent;
            const Particle* grandparent;
            std::vector<Particle> daughters;
            std::vector<Particle*> daughtersPtr;
            std::unique_ptr<angular::AngularCorrelation> angularCorrelation;
            std::unique_ptr<WeightCalculator> weightCalculator;
            std::unique_ptr<FinalStateGenerator> fsGenerator;            
            TVector3 Ox, Oy, Oz; // orthogonal vectors which define the Nuclear Process Coordinate System (NPCS)
            int orbitalL;
            Z_AXIS coordinateSystem;
            double smallNumber;

            void determineOxOyOz();
            void scale( std::vector<TLorentzVector>& );
    };
}

#endif	/* NBODYDECAY_H */

