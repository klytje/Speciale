#ifndef SIMX_BEAM_H
#define SIMX_BEAM_H

#include <TVector3.h>
#include <memory>
#include <TF1.h>

namespace ROOT {namespace Math {class DistSampler;}}

namespace simX {
    class Particle;

    class Beam {

        public:
            using FunctionPtr = std::unique_ptr<TF1>;
            Beam( double energy0, double x0, double y0, double z0, double theta, double phi,
                  FunctionPtr fEn=nullptr, FunctionPtr fXY=nullptr, FunctionPtr fTheta=nullptr, FunctionPtr fPhi=nullptr, bool perNucleon=false);

            Beam(Beam&&);
            ~Beam();

            /**
             * Returns energy the energy sampled from the user specified distribution.
             * @param A Nucleon number. If the beam is set to be per nucleon, the beam energy is multiplied by A before
             *        it is returned
             */
            double sampleEnergy(int A);

            /**
             * Return the nominal energy ie. without sampling the distribution.
             * @param A Nucleon number. If the beam is set to be per nucleon, the energy is multiplied by A before
             *        it is returned
             */
            double getNominalEnergy(int A);

            /**
             * If the beam is per nucleon, energies are will be multiplied by the
             * specified nucleon number A (parameter is relavant functions).
             * @return Is beam energy per nucleon
             */
            bool isPerNucleon() {
                return perNucleon;
            }

            /**
             * Sets the nominal energy. Notice, if the beam is per nucleon, this energy should also be per nucleon.
             * @param E
             */
            void setNominalEnergy(double E);

            /**
            * Returns XY position of beam particle
            * @return X and Y sampled from the user-specified distribution.
            */
            void sampleXY( double& x, double& y);

            /**
             * Returns the nominal XY position.
             */
            void nominalXY(double& x, double& y);

            /**
            * Returns direction of beam.
            * @return The direction is sampled from the user-specified angular distribution.
            */
            TVector3 sampleDirection();

            const TVector3& nominalDirection();

            double getZ0() {return z0;}

        private:
            double energy0, x0, y0, z0, theta, phi;
            bool samplerEnIsInitialized, samplerXYIsInitialized, samplerThetaIsInitialized, samplerPhiIsInitialized, perNucleon;
            FunctionPtr fEn, fXY, fTheta, fPhi;
            std::unique_ptr<ROOT::Math::DistSampler> samplerEn, samplerXY, samplerTheta, samplerPhi;
            TVector3 vNominalDirection;
    };
}

#endif	/* BEAM_H */
