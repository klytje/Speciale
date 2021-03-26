
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <simX/Beam.h>
#include <ausa/util/memory>
#include <ausa/constants/Mass.h>
#include <simX/Target.h>
#include <simX/EventSimulator.h>
#include <simX/propagator/NonIonizingPropagator.h>
#include <simX/PhysicsContainer.h>
#include <ausa/geometry/Box.h>
#include "simX/Particle.h"
#include "simX/parser/ReactionParser.h"


using namespace simX;
using namespace std;
using namespace simX::parser;
using namespace AUSA::EnergyLoss;

SUITE(PhysicsContainerTest) {

    Target buildTarget(int n) {
        std::vector<Layer> v;
        auto silicon = Material::predefined("Silicon");
        for (int i = 1; i <= n; i++) {
            v.emplace_back(silicon, make_unique<AUSA::Geometry::Box>(100.,100.,100.,TVector3(2.,-3,(i-1)*100.),TVector3(0,0,-1),TVector3(0,1,0)), i==n);
        }

        return {move(v), TVector3(2.,-3,100.)};
    }
    
    void setProp(ProcessChain& c) {
        auto prop = std::make_shared<propagator::NonIonizingPropagator>();
        c.getBeam().setPropagator(prop);

        for (auto& i : c.getTree()) {
            for (auto& j : i->getDaughters()) {
                j->setPropagator(prop);
            }
        }
    }

    class SetupFixture {
        public:
            SetupFixture()
                : beam(1E3, 0, 0, -10, 0, 0), target(buildTarget(1))
            {
                prop = std::make_shared<propagator::NonIonizingPropagator>();
            }

            ReactionParser parser;
            Beam beam;
            Target target;
            std::shared_ptr<propagator::NonIonizingPropagator> prop;
    };

    TEST_FIXTURE(SetupFixture, CreatePhysicsContainerAndWriteToDisk) {
        // simple process chain
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                He4
                Be8
            }
            )");
        // create physics container
        PhysicsContainer pCont("tmp.root");
    }
    

    TEST_FIXTURE(SetupFixture, FillPhysicsContainer) {
        // simple process chain
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                He4
                Be8
            }
            )");
        // set propagators
        setProp(*c);
        // create physics container
        PhysicsContainer* pCont = new PhysicsContainer("tmp.root");
        // create event simulator
        EventSimulator sim(beam, target, *c);
        // simulate 100 events and fill container
        for (int i=0; i<100; i++) {
            auto event = sim.run();
            pCont->fill(event);
        }
        // write to disk 
        delete pCont;

        // check that tree contains 100 events
        TFile *f = new TFile("tmp.root");
        TTree *tree = (TTree*)f->Get("phys");
        int N = tree->GetEntries();
        CHECK_EQUAL( 100, N );
    }

    TEST_FIXTURE(SetupFixture, DataWrittenToDiskIsCorrect) {
        // simple process chain
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                He4
                Be8
            }
            )");
        // set propagators
        setProp(*c);
        // create physics container
        PhysicsContainer* pCont = new PhysicsContainer("tmp.root");
        // create event simulator
        EventSimulator sim(beam, target, *c);
        // simulate 1 event
        auto event = sim.run();

        // Store energy for comparison
        vector<float> E;
        for (auto& p : event) E.push_back(p->getKineticEnergyLab());

        pCont->fill(event);

        // write to disk
        delete pCont;

        // check that tree contains 100 events
        TFile *f = new TFile("tmp.root");
        TTree *tree = (TTree*)f->Get("phys");

        vector<float> fileE(E.size());
        tree->SetBranchAddress("ENERGY", fileE.data());

        tree->GetEntry(0);

        CHECK_ARRAY_CLOSE(E, fileE, E.size(), 1E-6);
    }

    TEST_FIXTURE(SetupFixture, TFileWrapperCanBeUsedAsSourceForPhysicsContainer) {
        // simple process chain
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                He4
                Be8
            }
            )");
        // set propagators
        setProp(*c);
        // create physics container
        auto wrapper = std::make_shared<AUSA::TFileWrapper>("tmp.root", "RECREATE");
        PhysicsContainer* pCont = new PhysicsContainer(wrapper);
        // create event simulator
        EventSimulator sim(beam, target, *c);
        // simulate 100 events and fill container
        for (int i=0; i<100; i++) {
            auto event = sim.run();
            pCont->fill(event);
        }
        // write to disk
        delete pCont;
        wrapper->close();

        // check that tree contains 100 events
        TFile *f = new TFile("tmp.root");
        TTree *tree = (TTree*)f->Get("phys");
        int N = tree->GetEntries();
        CHECK_EQUAL( 100, N );
        delete f;
    }
}
