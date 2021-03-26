
#include "simX/propagator/MCStragglingPropagator.h"
#include "simX/propagator/Util.h"
#include "simX/Particle.h"
#include "simX/Random.h"

#include "ausa/constants/Mass.h"
#include "ausa/constants/Constants.h"
#include "ausa/eloss/EnergyLossIntegrator.h"
#include "ausa/eloss/SRIM13Tabulation.h"
#include "ausa/eloss/Material.h"
#include "ausa/eloss/StoppingPowerInterpolator.h"
#include <ausa/eloss/Ion.h>
#include <ausa/util/memory>

#include <cassert>
#include <cmath>

using namespace std;
using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;


MCStragglingPropagator::MCStragglingPropagator(LossCalcFactory lossFactory, RangeCalcFactory rangeFactory)
 : lowECutOff(10),
   csMultiplier(0.04),
   energyStraggling(true),
   lossFactory(lossFactory),
   rangeFactory(rangeFactory)
{
    double zmin = 0.001;
    double mp = AUSA::Constants::isotopeMass(1,1);
    double zmax = 1E5 / getEnergyUnitKeV(1, mp, 1, mp); // 100 MeV protons in hydrogen
    Jz = make_unique<ScatteringIntegral>(zmin, zmax);

    // histogram
    hn = new TH1F("hn", "Relative number densities", 10, 0, 10);
    hn->SetDirectory(nullptr);
}


void MCStragglingPropagator::setEnergyStraggling(bool status)
{
    energyStraggling = status;
}


void MCStragglingPropagator::setCrossSectionMultiplier(double multiplier)
{
    csMultiplier = multiplier;
}


void MCStragglingPropagator::setLowEnergyCutOff(double cutOff)
{
    lowECutOff = cutOff;
}


void MCStragglingPropagator::propagate(const Layer& layer, Particle& part, double range, double * nioni) 
{
    if (range == -1) range = findMaxRange(layer, part);

    // multiplicative factor
    double C = sqrt(csMultiplier) * 0.5;

    // Z1, M1
    double Z1 = part.getZ();
    double M1 = part.getMass();

    // get material
    auto& mat = layer.getMaterial();
    
    // number densities of each target element (in mm^-3)
    std::vector<double> numberDensities;
    double ntot = 0; // total number density (in mm^-3)
    for (size_t i=0; i<mat.getStoichiometries().size(); i++) {
        double n = mat.getStoichiometries()[i] * 1E-3 * mat.getDensity() / mat.getMassInGrams();
        numberDensities.push_back(n);
        ntot += n;
    }

    // Make histogram of relative number densities
    size_t i = 0;
    hn -> Reset();
    for (auto& n : numberDensities) {
        hn -> Fill(i+0.5, n/ntot);
        i++;
    }

    // If we have not worked with the Layer (ie. Material) before create LossCalc
    // for the ion and for all the target elements, and cache them for the next time.
    
    // LossCalc for ion
    auto key = make_pair(&layer, &part);
    if (!lossMap.count(key)) lossMap[key] = lossFactory(layer, part);
    auto& calc_ion = lossMap[key];

    // RangeCalc for ion
    if (!rangeMap.count(key)) rangeMap[key] = rangeFactory(layer, part);
    auto& range_ion = rangeMap[key];

    // LossCalc for target nuclei
    recoilKeys.clear();
    for (size_t i=0; i<mat.getElements().size(); i++) {		
        int z = mat.getElements()[i];
		int a = (int)(mat.getAtomicWeights()[i] / AUSA::Constants::ATOMIC_MASS_UNIT); // average A
		Particle recoil(Ion(z, a));
		auto key = make_pair(&layer, &recoil);
        if (!lossMap.count(key)) lossMap[key] = lossFactory(layer, recoil);
        recoilKeys.push_back(key);
    }

    
    // low-energy cut-off
    double emin = max(lowECutOff, calc_ion->getCutOff());
        
    // start position
    auto pos_i = part.getPosition();

    // start direction
    auto dir_i = part.getDirectionLab();
    
    // start energy
    auto energy_i = part.getKineticEnergyLab();
    
    // mean free path
    meanFreePath = 1. / (TMath::Pi() * pow(C, 2.) * cbrt(ntot));

    // ionizing and non-ionizing energy loss
    ionizingEnergyLoss = 0;
    nonIonizingEnergyLoss = 0;
    
    // number of collisions
    collisionCounter = 0;

    // begin scattering process
    auto pos     = pos_i;
    auto dir     = dir_i;
    auto energy  = energy_i;
    bool inside  = true;
    double projD = 0;
    while (inside && energy > emin && projD <= range) 
    {    
		// Generate randum number t : 0-1
		double t1 = rnd();
		
		// Travelling distance
		double x = -std::log(1. - t1);
		double D = x * meanFreePath;

		// Update position
		pos += D * dir;
		
		// Distance travelled projected onto initial direction
		projD = std::abs((pos-pos_i).Dot(dir_i));
		
		// check that particle is still inside the layer
		inside = layer.isInside(pos, VOLUME_TOLERANCE);
		
		// if particle has exited the layer, back trace trajectory to point of exit
		if (!inside) {
		    // find point of intersection
            TVector3 intersect;
            double dist, thick;
		    layer.getIntersection(pos, -dir, intersect, thick, dist);
		    // back trace
            pos = intersect;
            D -= dist;
		}

		// Electronic energy loss
		double dEe = calc_ion -> getElectronicEnergyLoss(energy, D);
		
		// Substract electronic energy loss
		energy -= dEe;
		
	    // increment ionizing energy loss
		ionizingEnergyLoss += dEe;
		
		// if ion was brought to rest, correct path length and terminate loop
		if (energy <= emin) {
		    double r = min(D, 1E-6*range_ion->getRange(energy + dEe));
		    pos += (r-D) * dir;
		    break;
		}

        // Electronic energy-loss straggling
        double eav = energy + dEe/2;
        double std = getEnergyStdDev(layer, part, D, eav);		

        // sample gaussian
        double straggl = 0;
        if (energy > 0. && energyStraggling) { 
            straggl = std * gaus(); 
            while (energy - straggl < 0) straggl = std * gaus();
        }
        
        // add straggling
        energy             -= straggl;
        dEe                += straggl;
        ionizingEnergyLoss += straggl;
		
		// if ion has exited the material, terminate loop
		if (!inside) break;
		
		// Determine scattering centre
		int elemId = hn -> GetRandom();

        // Z2, M2, and number density
		double Z2  = mat.getElements()[elemId];
		double M2  = mat.getAtomicWeights()[elemId];
		double n   = numberDensities[elemId];

        // some derived quantities
		double a    = getScreeningRadius(Z1, Z2);
		double r0   = C / cbrt(n); 
		double Jtot = pow(r0/a, 2.);
		
		// Reduced energy
		double eUnit = getEnergyUnitKeV(Z1, M1, Z2, M2);
		double epsilon = energy / eUnit;
		
		// Generate randum number t : 0-1
		double t2 = rnd();

		// Reduced c.m. scattering angle (Eq. 10)
		double Jeps = Jz -> eval(epsilon); 
		double j = Jeps + (t2 - 1.) * Jtot;
		double z = Jz -> invert(j, 0, epsilon);

		// Energy transfered to target atom in collision
		double mu = (M1 * M2) / (M1 + M2);
		double depsn = 4. * mu / (M1 + M2) * pow(z, 2.) / epsilon;
		double dEn = depsn * eUnit;
		
		// Subtract nuclear energy loss from ion energy
		energy -= dEn;

		// Calculate electronic and nuclear energy loss for recoiling nucleus
		// assuming that it comes to a full stop
        auto& calc_recoil = lossMap[recoilKeys[elemId]];
        double dEe_recoil = calc_recoil -> getElectronicEnergyLoss(dEn, 1E3); // thickness = 1 m
        double dEn_recoil = dEn - dEe_recoil;

	    // increment ionizing and non-ionizing energy loss
		ionizingEnergyLoss    += dEe_recoil;
		nonIonizingEnergyLoss += dEn_recoil;

		// Sample azimuthal lab angle (Eq. 12)
		double t3 = rnd();
		double psi = 2. * TMath::Pi() * t3;
		double cp = cos(psi);
		double sp = sin(psi);
		
		// Polar lab angle (Eq. 13)
		double ct = (1. - 2.*mu/M1 * pow(z/epsilon, 2.)) 
		            * 1./sqrt(1. - 4.*mu/(M1+M2) * pow(z/epsilon, 2.));
		double st = sqrt(1. - pow(ct, 2.));

		// Update direction of motion
		TVector3 v(st*cp, st*sp, ct);
		v.RotateUz(dir);
		dir = v.Unit();
		
		// increment counter
		collisionCounter++;
	}
	
	// if energy is below cut-off, use SRIM dE/dx tabulation to calculate dEe and dEn
	// for the last bit
	if (energy <= emin) {
	    double dEe = calc_ion -> getElectronicEnergyLoss(energy, 1E3);  // thickness = 1 m
	    double dEn = energy - dEe;
	    // increment ionizing and non-ionizing energy loss
		ionizingEnergyLoss    += dEe;
		nonIonizingEnergyLoss += dEn;
	    energy = 0;
	}
    	
	// update particle properties
    part.setFourMomentumLab(energy, dir);  	
    part.setPosition(pos);
    
    // return non-ionizing energy loss
    if (nioni != nullptr) *nioni = nonIonizingEnergyLoss;
}


double MCStragglingPropagator::getIonizingEnergyLoss() 
{
    return ionizingEnergyLoss;
}


double MCStragglingPropagator::getNonIonizingEnergyLoss() 
{
    return nonIonizingEnergyLoss;
}


size_t MCStragglingPropagator::getNumberOfCollisions()
{
    return collisionCounter;
}


double MCStragglingPropagator::getMeanFreePath()
{
    return meanFreePath;
}


double MCStragglingPropagator::getLowEnergyCutOff()
{
    return lowECutOff;
}


double MCStragglingPropagator::getCrossSectionMultiplier()
{
    return csMultiplier;
}


double MCStragglingPropagator::getEnergyUnitKeV(double Z1, double M1, double Z2, double M2) 
{
    double alpha = AUSA::FINE_STRUCTURE_CONSTANT;
    double me = AUSA::Constants::ELECTRON_MASS;

    double mu = M1 * M2 / (M1 + M2);

    // factor in Eq. (2) [keV^-1]
    double f = mu/M1 * 1./(Z1*Z2) * 1./pow(alpha, 2.) * 0.8853/me * 1./sqrt(cbrt2(Z1) + cbrt2(Z2)); 
    
    return 1./f;
}


double MCStragglingPropagator::getScreeningRadius(double Z1, double Z2)
{
    return  0.8853 * AUSA::BOHR_RADIUS * 1./sqrt(cbrt2(Z1) + cbrt2(Z2));
}
