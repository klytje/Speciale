//
// Created by munk on 26-08-15.
//

#include "simX/propagator/GaussianStragglingPropagator.h"
#include "simX/propagator/Util.h"
#include "simX/Particle.h"
#include "simX/Random.h"

#include "ausa/constants/Mass.h"

#include <cassert>

using namespace simX;
using namespace simX::propagator;
using namespace std;

GaussianStragglingPropagator::GaussianStragglingPropagator(
    std::shared_ptr<ParticlePropagator> inner, bool angularStraggling, bool energyStraggling)
 : screen(SCREENING::THOMAS_FERMI), 
   inner(inner),
   halfScattAngle(0),
   nsig(10), // how many sigmas we want to sample
   angularStraggling(angularStraggling),
   energyStraggling(energyStraggling)
          
{
    assert(inner != nullptr);
    
    // initialize universal small-angle scattering function
    double xmin = 0;
    double xmax = nsig / sqrt(2.);
    xmax = TMath::Min(xmax, TMath::Pi());
    scattFunc = TF1("scattFunc", "x*exp(-x*x)", xmin, xmax); 
    
    ws.resize(4, 0);
}


void GaussianStragglingPropagator::setScreeningModel(GaussianStragglingPropagator::SCREENING screening)
{
    screen = screening;
}


void GaussianStragglingPropagator::propagate(const Layer& layer, Particle& part, double range, double * nonIonizingEloss) {

    // reset half-scattering angle
    halfScattAngle = 0;    

    // Initial particle energy
    auto ei = part.getKineticEnergyLab();

    // Initial particle position
    auto pi = part.getPosition(); 
    
    // Let inner "normal" propagator do its stuff first
    inner -> propagate(layer, part, range);

    // Final particle energy
    double ef = part.getKineticEnergyLab();
    
    // if particles comes to a full stop, end calculation here (i.e. no straggling for fully stopped particles)
    if (ef == 0) return;
    
    // Final position
    auto pf = part.getPosition();

    // Distance propagated
    auto dist = (pf - pi).Mag();
    
    // Energy straggling
    if (energyStraggling) doEnergyStraggling(layer, part, dist, ei, ef);

    // Angular straggling
    if (angularStraggling) doAngularStraggling(layer, part, dist, ei, ef);

    // All energy loss is assumed to be ionizing
    if (nonIonizingEloss != nullptr) *nonIonizingEloss = 0;
}


double GaussianStragglingPropagator::getHalfScatteringAngle() const
{
    return halfScattAngle * 180. / TMath::Pi();
}


void GaussianStragglingPropagator::doAngularStraggling(const Layer& layer, Particle& part, double dist, double ei, double ef) 
{
    // Old direction
    auto dir = part.getDirectionLab();

    // Mean energy
    double eav = (ei + ef) / 2.;

    // Calculate half-scattering angle based on empirical parameterisations
    halfScattAngle = getHalfScatteringAngle(layer, part, dist, eav);

    // 1-sigma spread
    auto sigma = 2./2.35 * halfScattAngle;

    // Sample angles
    double xmin = 0;
    double xmax = 0;
    scattFunc.GetRange(xmin, xmax);
    double alpha = sqrt(2.) * sigma * scattFunc.GetRandom(xmin, xmax); // scattering angle
    double beta  = 2. * TMath::Pi() * rnd();                           // azimuthal angle

    // deflection vector
    TVector3 newDir(0,0,1);
    newDir.SetTheta(alpha);
    newDir.SetPhi(beta);

    // rotate deflection vector to obtain new direction
    // (not necessary if particle is moving along z axis)
    if (dir.Z() < 1) {
        TVector3 zAxis(0,0,1);
        TVector3 rotAxis = zAxis.Cross(dir);
        double rotAngle = dir.Theta();
        newDir.Rotate(rotAngle, rotAxis);    
    }
    
    // set new direction
    part.setFourMomentumLab(part.getKineticEnergyLab(), newDir);    
}


double GaussianStragglingPropagator::getHalfScatteringAngle(const Layer& layer, Particle& part, double dist, double energy)
{
    // particle properties
    double Z1 = part.getZ();
    double A1 = part.getA();
    double TA = energy / A1 / 1000.;  // kinetic energy per nucleon in MeV

    // get material
    auto& mat = layer.getMaterial();
    
    // loop over elements
    double taupSum = 0;
    double aapSum = 0;
    for (size_t i=0; i<mat.getNumberOfElements(); i++) 
    {
        // target properties
        double Z2 = mat.getElements()[i];
        double A2 = mat.getAtomicWeights()[i] / AUSA::Constants::ATOMIC_MASS_UNIT;

        //  double amuMeV = AUSA::CONSTANTS::ATOMIC_MASS_UNIT / 1000;
        //  double v = sqrt(1 - pow(1 + TA/amuMeV, -2.));  // speed in units of c (relativistic formula)
        //  double aB = Z1 * Z2 / (137 * v) ;
        //  cout << "Born parameter (nonrelativistic approximation): " << aB << endl;
        
        // Reduced kinetic energy of beam particle (epsilon_p)
        double Ep   = TA * A1;
        double Bpt  = 0.961 * Z1 * Z2 / cbrt(A2);
        double epsp = Ep / Bpt;

        // ratio of alpha-tilde to alpha for beam particle (Eq. 12)
        double Zpt = sqrt(cbrt2(Z1) + cbrt2(Z2));
        double aap = 15.63 * epsp / (cbrt(A2) * Zpt);

        // distance traveled in ug/cm2 (for this element only)
        double rx = (dist * 1E-1) * (mat.getDensity() * 1E6) * mat.getMassFractions()[i];

        // tau for beam particle
        double taup = 41.5 * rx / (A2 * pow(Zpt, 2));
        
        // sum
        taupSum += taup;

        // sum of aap ratios, weighted by taup
        aapSum += aap * taup;
    }

    // calculate alpha-tilde parameter        
    double a12Tilde = getAlphaTilde(taupSum);
    
    // average aap
    double aap = aapSum / taupSum;
    
    // half-scattering angle in radians
    return a12Tilde / aap / 1000.;
}


double GaussianStragglingPropagator::getAlphaTilde(double tau)
{
    double a12[4];

    for (size_t n=0; n<4; n++) 
    {
        double m=0, Cm=0;
        if (n == 0) { // tau<0.1
            if (screen == SCREENING::THOMAS_FERMI) {
                m  = 0.311;
                Cm = 1.05;
            }
            else if (screen == SCREENING::LENZ_JENSEN) {
                m  = 0.191;
                Cm = 3.45;
            }
        }
        else if (n == 1) { // Sigmund 1<tau<5
            m   = 0.5;
            Cm  = 0.25;
        }            
        else if (n == 2) { // Anne global
            m = 0.89;
            Cm = 0.92;
        }
        else if (n == 3) { // Anne tau>1e3
            m = 0.91;
            Cm = 1.00;
        }
        
        a12[n] = Cm * pow(tau, 1./(2.*m));
    }

    double res = 0;
    auto& f = getInterpolationFactors(tau);
    for (size_t n=0; n<4; n++) res += f[n] * a12[n];

    return res;
}


/* interpolation function to connect different tau-domains */
vector<double>& GaussianStragglingPropagator::getInterpolationFactors(double tau) {
    // reset ws
    for (int n=0; n<4; n++) ws[n] = 0;
    double d;
    double x = log10(tau);
    // different tau-domains
    if (tau<=0.1){ // <0.1
        ws[0] = 1;
    }
    if (tau>0.1 && tau<=1) { // 0.1-1
        d = 1.0;
        ws[0] = (log10(1)-x)/d;
        ws[1] = (x-log10(0.1))/d;
    }
    if (tau>1 && tau<=5) { // 1-5
        ws[1] = 1;
    }
    if (tau>5 && tau<=40) { // 5-40
        d = log10(40)-log10(5);
        ws[1] = (log10(40)-x)/d;
        ws[2] = (x-log10(5))/d;
    }
    if (tau>40 && tau<=500) { // 40-500
        ws[2] = 1;
    }
    if (tau>500 && tau<=1e3) { // 500-1e3
        d = log10(1e3)-log10(500);
        ws[2] = (log10(1e3)-x)/d;
        ws[3] = (x-log10(500))/d;
    }
    if (tau>1e3) { // 1e3-1e6
        ws[3] = 1;
    }
    return ws;
}


void GaussianStragglingPropagator::doEnergyStraggling(const Layer& layer, Particle& part, double dist, double ei, double ef) 
{
    // average energy
    double eav = (ei + ef) / 2.;

    // std
    double std = getEnergyStdDev(layer, part, dist, eav);
       
    // sample gaussian
    double e = -1; 
    while ( e < 0. ) e = ef + std * gaus(); 

    // update energy
    part.setFourMomentumLab(e, part.getDirectionLab());    
}


