//
// Created by munk on 06-08-15.
//

#include "simX/propagator/Util.h"
#include "simX/Layer.h"
#include "simX/Particle.h"
#include <ausa/constants/Mass.h>
#include <cmath>

namespace simX 
{
    namespace propagator 
    {
        double findMaxRange(const Layer& layer, const Particle& p) 
        {
            auto& pos = p.getPosition();
            auto& dir = p.getDirectionLab();

            TVector3 i;
            double d, r;
            layer.getIntersection(pos, dir, i, r, d);
            return r;
        }

        double getEnergyStdDev(const Layer& layer, Particle& part, double dist, double energy) 
        {
            // particle properties
            double Z1 = part.getZ();
            double A1 = part.getA();
            double TA = energy / A1 / 1000.;  // kinetic energy per nucleon in MeV

            // get material
            auto& mat = layer.getMaterial();

            // loop over elements
            double OmegaSumOfSquares = 0;
            for (size_t i=0; i<mat.getNumberOfElements(); i++) 
            {
                // target properties
                double Z2 = mat.getElements()[i];
                double A2 = mat.getAtomicWeights()[i] / AUSA::Constants::ATOMIC_MASS_UNIT;

                // distance traveled in ug/cm2 (for this element only)
                double rx = (dist * 1E-1) * (mat.getDensity() * 1E6) * mat.getMassFractions()[i];

                // Bohr straggling
                double OmegaB = sqrt(0.157 * rx * pow(Z1, 2) * Z2 / A2);
                
                // k factor
                double k = 1.1 + 0.47*std::log10(TA);

                // Straggling (keV)        
                double Omega = k * OmegaB;
                
                OmegaSumOfSquares += pow(Omega, 2);
            }

            return sqrt(OmegaSumOfSquares);
        }

    }
}
