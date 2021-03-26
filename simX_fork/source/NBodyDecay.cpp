
#include "simX/NBodyDecay.h"
#include "simX/angular/AngularCorrelation.h"
#include "simX/generator/NBodyPhaseSpace.h"

#include "ausa/util/memory"
#include <assert.h>
#include <simX/Logger.h>

using namespace std;
using namespace simX;
using namespace simX::angular;

namespace {
    /**
     * Small number
     */
    double EPS = 1E-6;
}

NBodyDecay::NBodyDecay( Particle& parent, vector<Particle> theDaughters, std::unique_ptr<FinalStateGenerator> generator, unique_ptr<AngularCorrelation> ac, unique_ptr<WeightCalculator> wCalc, int L, Z_AXIS axis )
 :  parent(parent),
    daughters(theDaughters),
    MULT(theDaughters.size()),
    angularCorrelation(move(ac)),
    weightCalculator(move(wCalc)),
    fsGenerator(move(generator)),
    orbitalL(L),
    coordinateSystem(axis)
{
    // set parent
//    for (int i=0; i<MULT; i++) daughters[i].setParent(parent);

    // set daughters
    for (int i=0; i<MULT; i++) daughtersPtr.push_back( &daughters[i] );

    // get pointer to grandparent
    grandparent = parent.getParent();

    // default NuclearProcess Coordinate System (NPCS)
    Ox = TVector3(1,0,0);
    Oy = TVector3(0,1,0);
    Oz = TVector3(0,0,1);
    
    // Difference in rest mass (not including excitation energies) between
    // initial and final states
    mDiff = parent.getMass();
    for (auto& d : daughters) mDiff -= d.getMass();
}

NBodyDecay::NBodyDecay(NBodyDecay&& other) :
    MULT(other.MULT),
    mDiff(other.mDiff),
    parent(other.parent),
    grandparent(other.grandparent),
    daughters(move(other.daughters)),
    daughtersPtr(move(other.daughtersPtr)),
    angularCorrelation(move(other.angularCorrelation)),
    weightCalculator(move(other.weightCalculator)),
    fsGenerator(move(other.fsGenerator)),
    Ox(other.Ox),
    Oy(other.Oy),
    Oz(other.Oz),
    orbitalL(other.orbitalL),
    coordinateSystem(other.coordinateSystem)
{
    // Nothing
}

NBodyDecay::~NBodyDecay() {
    // Do nothing
}

NBodyDecay NBodyDecay::withDefaultGenerator(Particle& parent, std::vector<Particle> theDaughters,
                                              unique_ptr<AngularCorrelation> ac,
                                              unique_ptr<WeightCalculator> wCalc) {
    NBodyDecay n(parent, theDaughters, nullptr, move(ac), move(wCalc), 0);
    n.setFinalStateGenerator(std::make_unique<NBodyPhaseSpace>(parent, n.getDaughters()));
    return n;
}

void NBodyDecay::runProcess() {
    static auto log = log::getLogger("NBodyDecay");

    assert(fsGenerator!=nullptr);
    assert(angularCorrelation!=nullptr);

    // Place daughters @ parent position
    for (auto& d : daughters)
        d.setPosition(parent.getPosition());

    // Determine the 3 axes of the coordinates system
    determineOxOyOz();

    // Generate final state
    auto& PD = fsGenerator->getFourMomenta();
    
    // If needed, scale to reproduce current excitation energy of parent
    // (This is cheating, but usually a decent approximation if the 
    // the excitation energy does not vary too much compared to the
    // the total kinetic energy of the daugthers)
    if (parent.hasFiniteWidth()) scale(PD);
    
    // Determine boost vector from parent's rest frame -> lab
    auto b = parent.getFourMomentumLab().BoostVector();

    // Sample rotation angles
    double theta=0, phi=0;
    angularCorrelation->sampleAngles(theta, phi);
    log->debug("Rotate theta: {} phi: {}", theta, phi);

    // Loop over daughters
    for (int i=0; i<MULT; i++) {
        // Rotate TVector3 component
        TLorentzVector& pI = PD[i];



        pI.RotateY(theta);
        pI.RotateZ(phi);
        log->debug("Daughter {} four momentum (p,E) = ({}, {}, {}, {})", daughters[i].getName(), pI.E(), pI.X(), pI.Y(), pI.Z());

        // Change coordinate system (x,y,z)->(Ox,Oy,Oz)
        TVector3 p = pI.Px()*Ox + pI.Py()*Oy + pI.Pz()*Oz;
        pI.SetVect(p);
    
        // Boost to LAB and update
        pI.Boost(b);
	    daughters[i].setFourMomentumLab(pI);
    }
}

vector<Particle*>&NBodyDecay::getDaughters() {
    return daughtersPtr;
}

const std::vector<Particle *>& NBodyDecay::getDaughters() const {
    return daughtersPtr;
}

void NBodyDecay::determineOxOyOz() {
    // Oz
    Oz.SetXYZ(0,0,1);
    TVector3 a = parent.getDirectionRFOP();
    if (a.Mag() >= EPS) Oz = a;
    // Oy
    TVector3 b(0,0,1);
    if (grandparent != nullptr) b = grandparent -> getDirectionRFOP();
    TVector3 c = b.Cross(Oz);
    if (c.Mag() >= EPS) Oy = c;
    else {
        c = TVector3(-1,0,0).Cross(Oz);
        if (c.Mag() < EPS) c = TVector3(0, 1, 0);
    }
    Oy = c.Unit();
    // Ox
    Ox = Oy.Cross(Oz).Unit();
    // Switch to Alan Shotter's coordinate system, if requested
    if (coordinateSystem==Z_AXIS::ORTHOGONAL) {
        TVector3 temp = Oz;
        Oz = -Oy;
        Oy = temp;    
    }
}


void NBodyDecay::setFinalStateGenerator(std::unique_ptr<FinalStateGenerator> generator) {
    fsGenerator = move(generator);
}

void NBodyDecay::scale( vector<TLorentzVector>& pd ) {
    // initial energy
    double ei = parent.getMass() + parent.getExcitationEnergy();
    // final energy
    double ef = 0;
    // final kinetic energy
    double kf = 0;
    for (auto& p : pd) {
        ef += p.E();
        kf += p.E() - p.Mag();
    }
    // scaling factor to be applied to kinetic energies
    double sk = 1;
    if (ef>0) sk += (ei-ef)/kf;
    // scaling factor to be applied to momenta
    // (guarantees that m^2=E^2-p^2 is fulfilled, but only in the
    // non-relativistic limit)
    double sp = sqrt(sk);    
    // scale kinetic energies and momenta
    for (auto& p : pd) {
        double m = p.Mag();
        double k = p.E() - m;
        k *= sk;
        p.SetE(m+k);
        TVector3 p3 = p.Vect();
        p3 *= sp;
        p.SetVect(p3);
    }
}


void NBodyDecay::setWeightCalculator(unique_ptr<WeightCalculator> wCalc) {
    weightCalculator = move(wCalc); 
}


double NBodyDecay::getWeight() const {
//    assert(weightCalculator!=nullptr);
    double w = 1;
    if (weightCalculator!=nullptr) w = weightCalculator -> getWeight();
    return w;
}


int NBodyDecay::getNumberOfNucleons() const {
    return parent.getA();
}

int NBodyDecay::getNumberOfProtonsInInitialState() const {
    return parent.getZ();
}


