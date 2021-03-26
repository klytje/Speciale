
#include "simX/EventSimulator.h"
#include "simX/propagator/Vacuum.h"
#include "simX/Random.h"
#include "ausa/util/stream.h"
#include "simX/Logger.h"

#include <TVector3.h>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <simX/CompoundFormation.h>

using namespace std;
using namespace simX;

namespace {
    enum {
        NOMIAL_BEAM = 0,
        SAMPLE_BEAM = 1
    };
}


EventSimulator::EventSimulator( Beam& beam, Target& target, ProcessChain& c )
 :  beam(beam), 
    beamPart(c.getBeam()),
    chain(c),
    event(-1)
{
    assert(!target.getLayers().empty());

    // Identify final-state particles
    finalStates = chain.findFinalStates();
    
    // Order target layers according to increasing distance as beam sees it.
    // (Must be done before calling any of the propagation routines!)
    orderTargetLayers(target);

    // Propagate beam (with parameters set to average values) 
    // halfway into the active layer of the target. Then determine 
    // the nominal excitation energy of the compound nucleus.
    // Only do this for process chains that begin with compound formation.
    auto first = chain.getFirstProcess();
    if (dynamic_cast<CompoundFormation*>(first)) {
    
        // initialize beam with nominal values
        sampleBeam(NOMIAL_BEAM);
        
        // propagate halfway into active layer
        propagateBeam(0.5);
        
        // simulate compound formation
        first -> runProcess();
        
        // set nominal excitation energy 
        double ex = first -> getDaughters()[0] -> getExcitationEnergy();
        first -> getDaughters()[0] -> setNominalExcitationEnergy( ex );
    }
        
    // Create excitation-energy sampler
    // (Important that this happens after the nominal excitation energy of
    // the compound nucleus has been determined!)
    exSampler = unique_ptr<ExcitationSampler>( new ExcitationSampler( chain ) );

    {
        auto log = log::getLogger("EventSimulator");
        double x, y;
        auto& dir = beam.nominalDirection();
        beam.nominalXY(x,y);
        log->debug("Beam nominal position: ({}, {}, {})", x, y, beam.getZ0());
        log->debug("Beam nominal direction: ({}, {}, {})", dir.x(), dir.y(), dir.z());
    }
}


EventSimulator::~EventSimulator() 
{}


PhysicsEvent EventSimulator::run() {
    static auto log = log::getLogger("EventSimulator");
    log->debug("Simulating event {}", ++event);

    // sample excitation energies of intermediate resonances
    exSampler -> sampleExcitationEnergies();

    // initialize the beam particle
    sampleBeam(SAMPLE_BEAM);
    
    // propagate the beam to the reaction/decay site
    auto react = propagateBeam(simX::rnd());

    // If we did not hit the active layer then we stop here
    if (!react) return PhysicsEvent{empty};

    // run all processes
    for (NuclearProcess & d : chain) {
        d.runProcess();
    }

    // propagate out through layers
    for (auto p : finalStates) {
        propagateParticleOut(*p);
    }

    return PhysicsEvent(finalStates);
}


void EventSimulator::sampleBeam( int mode ) {
    double x, y, T;
    TVector3 dir;
    if (mode == NOMIAL_BEAM) {
        beam.nominalXY(x,y);
        dir = beam.nominalDirection();
        T = beam.getNominalEnergy(beamPart.getA());
    }
    else {
        beam.sampleXY(x, y);
        dir = beam.sampleDirection();
        T = beam.sampleEnergy(beamPart.getA());
    }
    beamPart.setFourMomentumLab(T,dir);
    beamPart.setPosition(TVector3(x,y,beam.getZ0()));
}

bool EventSimulator::propagateBeam( double f ) {
    static auto log = log::getLogger("EventSimulator");
    // Propagate from -10cm to active layer
    auto i = propagator::propagateInVacuum(layers, beamPart);
    if (i == -1) return false;

    auto reaction = false;
    log->debug("Beam energy before target: {}", beamPart.getKineticEnergyLab());
    for (size_t j = 0; j <= indexActiveLayer; ++j) {
        auto layer = layers[j];
        TVector3 intersection;
        double t, d;
        bool hit = layer->getIntersection(beamPart.getPosition(), beamPart.getDirectionLab(), intersection, t, d);
        if (!hit) return reaction;

        t*= layer->isActive() ? f : 1;
        reaction |= layer->isActive();

        beamPart.propagate(*layer, beamPart, t);
        log->debug("Beam energy after layer {}: {}", layer->getMaterial().getName(), beamPart.getKineticEnergyLab());
    }
    auto& p = beamPart.getPosition();
    log->debug("Reaction site: ({}, {}, {})", p.x(), p.y(), p.z());
    return reaction;
}

// Purpose: Create a vector with layers ordered according distance as beam sees it.
void EventSimulator::orderTargetLayers(const Target& target) {
    using Pair = std::pair<double,const Layer*>;

    sampleBeam(NOMIAL_BEAM); // Place beam far away

    // Vector with ordered pair of dist/layer
    vector<Pair> v;
    v.reserve(target.getLayers().size());

    // Tmp vars for intersection algo
    TVector3 intersection;
    double t, dist;

    // Compute distance to each.
    for (auto& layer : target.getLayers()) {
        auto hit = layer.getIntersection(beamPart.getPosition(), beamPart.getDirectionLab(), intersection, t, dist);
        v.emplace_back(dist, &layer);

        if (!hit) throw std::runtime_error("The nominal beam cannot hit target layer '" + layer.getMaterial().getName() + "'");
    }

    // Sort according to dist
    std::sort(begin(v), end(v), [](const Pair& p1, const Pair& p2) {
        return p1.first < p2.first;
    });

    assert(v[0].first >= 0 && "Beam cannot hit one or more of your target layers");

    // Replace. Remember which layer is active.
    indexActiveLayer = static_cast<size_t>(-1);
    for (size_t i = 0; i < v.size(); ++i) {
        auto l = v[i].second;
        layers.push_back(l);
        if (l->isActive()) {
            if (indexActiveLayer != -1) throw runtime_error("There must be only one active layer!");
            indexActiveLayer = i;
        }
    }
    if (indexActiveLayer == -1) throw runtime_error("There must be a active layer!");
}

void EventSimulator::propagateParticleOut(Particle &p) {
    static auto log = log::getLogger("EventSimulator");
    auto next = findNextLayer(p, nullptr);
    log->debug("Daughter '{}' energy at reaction site: {}", p.getName(), p.getKineticEnergyLab());
    while (next != nullptr && p.getKineticEnergyLab() > 0) {
        p.propagate(*next, p);
        log->debug("Daughter '{}' energy after {}: {}", p.getName(), next->getMaterial().getName(), p.getKineticEnergyLab());
        next = findNextLayer(p, next);
    }
}

const Layer *EventSimulator::findNextLayer(Particle &p, const Layer* old) {
    for (auto layer : layers) {
        if (layer != old && layer->isInside(p.getPosition(), VOLUME_TOLERANCE)) {
            return layer;
        }
    }
    return nullptr;
}
