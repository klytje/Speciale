
#include "simX/Beam.h"
#include "simX/Random.h"

#include <TMath.h>
#include <cmath>

#include "Math/DistSampler.h"
#include "Math/DistSamplerOptions.h"
#include "Math/Factory.h"

using namespace std;
using namespace simX;
using ROOT::Math::Factory;
using DistPtr = std::unique_ptr<ROOT::Math::DistSampler>;

namespace {
    const auto SAMPLER_TYPE = "Foam";

    DistPtr build1DSampler(const Beam::FunctionPtr& f) {
        if (!f) return DistPtr{};
        auto sampler = DistPtr(Factory::CreateDistSampler(SAMPLER_TYPE));
        sampler -> SetFunction( *f, 1 );
        double minEn;
        double maxEn;
        f -> GetRange(minEn, maxEn);
        sampler -> SetRange(minEn,maxEn);
        return sampler;
    }
}


Beam::Beam(Beam && other) :
energy0(other.energy0),
x0(other.x0),
y0(other.y0),
z0(other.z0),
theta(other.theta),
phi(other.phi),
samplerEnIsInitialized(other.samplerEnIsInitialized),
samplerXYIsInitialized(other.samplerXYIsInitialized),
samplerThetaIsInitialized(other.samplerThetaIsInitialized),
samplerPhiIsInitialized(other.samplerPhiIsInitialized),
samplerEn(move(other.samplerEn)),
samplerPhi(move(other.samplerPhi)),
samplerTheta(move(other.samplerTheta)),
samplerXY(move(other.samplerXY)),
vNominalDirection(other.vNominalDirection)
{

}

Beam::Beam( double energy0, double x0, double y0, double z0, double theta, double phi,
            FunctionPtr fEn, FunctionPtr fXY, FunctionPtr fTheta, FunctionPtr fPhi, bool perNucleon)
 :  energy0(energy0),
    x0(x0),
    y0(y0),
    z0(z0),
    theta(theta),
    phi(phi),
    samplerEnIsInitialized(false),
    samplerXYIsInitialized(false),
    samplerThetaIsInitialized(false),
    samplerPhiIsInitialized(false),
    perNucleon(perNucleon)
{
    samplerEn = build1DSampler(fEn);
    samplerTheta = build1DSampler(fTheta);
    samplerPhi = build1DSampler(fPhi);

    // XY sampler (polar coordinates!)
    if (fXY!=nullptr) {
        samplerXY = DistPtr(Factory::CreateDistSampler(SAMPLER_TYPE));
        samplerXY -> SetFunction( *fXY, 2 ); 
        double minXY[2], maxXY[2];
        fXY -> GetRange(minXY[0], maxXY[0]);
        minXY[1] = 0.; 
        maxXY[1] = 2.*TMath::Pi();
        samplerXY -> SetRange(minXY,maxXY);
    }

    Beam::fEn = move(fEn);
    Beam::fXY = move(fXY);
    Beam::fTheta = move(fTheta);
    Beam::fPhi = move(fPhi);

    vNominalDirection = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
}

Beam::~Beam() 
{/*We have forward declared TF1 + DistSampler, so this needs to be here.*/}

double Beam::sampleEnergy(int A) {
    double v = 0;
    if (samplerEn) {
        if (!samplerEnIsInitialized) {  
            bool ret = samplerEn -> Init();	
            if (!ret)
                throw invalid_argument("Error --- beam energy sampler could not be initialized");
            samplerEnIsInitialized = true;
        }
        samplerEn -> Sample(&v);
    }

    return (perNucleon ? (energy0 + v) * A: energy0 + v);
}

void Beam::nominalXY(double &x, double &y) {
    x = x0;
    y = y0;
}

void Beam::sampleXY( double& x, double& y) {
    x = x0;
    y = y0;

    if (samplerXY) {
        double v[2] = {0,0};
        if (!samplerXYIsInitialized) {
            bool ret = samplerXY -> Init();	
            if (!ret)
                throw invalid_argument("Error --- beam position sampler could not be initialized");
            samplerXYIsInitialized = true;
        }
        samplerXY -> Sample(v);

        double r = v[0];
        double q = v[1];
        x += r*cos(q);
        y += r*sin(q);
    }
}

TVector3 Beam::sampleDirection() {
    if (!samplerPhi && !samplerTheta) return vNominalDirection;

    double dTheta = 0;
    double dPhi = 0;
    if (samplerTheta) {
        if (!samplerThetaIsInitialized) {
            bool ret = samplerTheta -> Init();
            if (!ret)
                throw invalid_argument("Error --- beam angle sampler could not be initialized");
            samplerThetaIsInitialized = true;
        }
        samplerTheta -> Sample(&dTheta);
    }

    if (samplerPhi) {
        if (!samplerPhiIsInitialized) {
            bool ret = samplerPhi -> Init();
            if (!ret)
                throw invalid_argument("Error --- beam angle sampler could not be initialized");
            samplerPhiIsInitialized = true;
        }
        samplerPhi -> Sample(&dPhi);
    }
    double theta = Beam::theta + dTheta;
    double phi = Beam::phi + dPhi;
    return TVector3( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
}


const TVector3& Beam::nominalDirection() {
    return vNominalDirection;
}

double Beam::getNominalEnergy(int A) {
    return perNucleon ? energy0*A : energy0;
}

void Beam::setNominalEnergy(double E) {
    energy0 = E;
}



