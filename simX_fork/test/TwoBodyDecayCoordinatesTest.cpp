
#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include "ausa/util/memory"
#include <simX/angular/IsotropicAngularCorrelation.h>
#include <simX/angular/FixedAngularCorrelation.h>

#include "simX/Particle.h"
#include "simX/TwoBodyDecay.h"
#include "simX/angular/AngularCorrelationTF1.h"
#include "simX/parser/ReactionParser.h"
#include "simX/ProcessChain.h"

using namespace simX;
using namespace std;
using namespace simX::parser;
using namespace simX::angular;
using namespace AUSA::EnergyLoss;

SUITE(TwoBodyDecayCoordinatesTest) {

    class SetupFixture {
        public:
            SetupFixture() {

            }
        ReactionParser px;      
    };

    TEST_FIXTURE(SetupFixture, OzIsAsExpectedForParentAtRest) {
        Particle Lu151( Ion(71,151), 0.0 );
        TwoBodyDecay p( Lu151, Ion(1,1), Ion(70,150), std::make_unique<FixedAngularCorrelation>(0,0) );
        p.runProcess();
        auto& d = p.getDaughters();
        TVector3 dir1 = d[0] -> getDirectionLab();
        CHECK_EQUAL( 0., dir1.X() );
        CHECK_EQUAL( 0., dir1.Y() );
        CHECK_EQUAL( 1., dir1.Z() );
        TVector3 dir2 = d[1] -> getDirectionLab();
        CHECK_EQUAL( 0., dir2.X() );
        CHECK_EQUAL( 0., dir2.Y() );
        CHECK_EQUAL( -1., dir2.Z() );        
    }

    TEST_FIXTURE(SetupFixture, OyIsAsExpectedForParentAtRest) {
        Particle Lu151( Ion(71,151), 0.0 );
        TwoBodyDecay p( Lu151, Ion(1,1), Ion(70,150), std::make_unique<FixedAngularCorrelation>(90,90) );
        p.runProcess();
        auto& d = p.getDaughters();
        TVector3 dir1 = d[0] -> getDirectionLab();
        CHECK_CLOSE( 0., dir1.X(), 1E-9 );
        CHECK_CLOSE( 1., dir1.Y(), 1E-9 );
        CHECK_CLOSE( 0., dir1.Z(), 1E-9 );
        TVector3 dir2 = d[1] -> getDirectionLab();
        CHECK_CLOSE( 0., dir2.X(), 1E-9 );
        CHECK_CLOSE( -1., dir2.Y(), 1E-9 );
        CHECK_CLOSE( 0., dir2.Z(), 1E-9 );        
    } 

    TEST_FIXTURE(SetupFixture, OzIsAsExpectedForParentInXYPlane) {
        Particle Lu151( Ion(71,151), 0.0 );
        double ekin=1.e-6;
        TVector3 dir(1,1,0);
        Lu151.setFourMomentumLab(ekin,dir);
        TwoBodyDecay p( Lu151, Ion(1,1), Ion(70,150), std::make_unique<FixedAngularCorrelation>(0,0) );
        p.runProcess();
        auto& d = p.getDaughters();
        TVector3 dir1 = d[0] -> getDirectionLab();
        CHECK_CLOSE( 1./std::sqrt(2.), dir1.X(), 1E-3 );
        CHECK_CLOSE( 1./std::sqrt(2.), dir1.Y(), 1E-3 );
        CHECK_CLOSE( 0., dir1.Z(), 1E-3 );
        TVector3 dir2 = d[1] -> getDirectionLab();
        CHECK_CLOSE( -1./std::sqrt(2.), dir2.X(), 1E-3 );
        CHECK_CLOSE( -1./std::sqrt(2.), dir2.Y(), 1E-3 );
        CHECK_CLOSE( 0., dir2.Z(), 1E-3 );        
    }

    TEST_FIXTURE(SetupFixture, OzIsAsExpectedForParentAlongZAxis) {
        Particle Lu151( Ion(71,151), 0.0 );
        double ekin=1.e-6;
        TVector3 dir(0,0,1);
        Lu151.setFourMomentumLab(ekin,dir);
        TwoBodyDecay p( Lu151, Ion(1,1), Ion(70,150), std::make_unique<FixedAngularCorrelation>(0,0) );
        p.runProcess();
        auto& d = p.getDaughters();
        TVector3 dir1 = d[0] -> getDirectionLab();
        CHECK_CLOSE( 0., dir1.X(), 1E-9 );
        CHECK_CLOSE( 0., dir1.Y(), 1E-9 );
        CHECK_CLOSE( 1., dir1.Z(), 1E-9 );
        TVector3 dir2 = d[1] -> getDirectionLab();
        CHECK_CLOSE( 0., dir2.X(), 1E-9 );
        CHECK_CLOSE( 0., dir2.Y(), 1E-9 );
        CHECK_CLOSE( -1., dir2.Z(), 1E-9 );        
    }

    TEST_FIXTURE(SetupFixture, OxIsAsExpectedForParentAlongZAxis) {
        Particle Lu151( Ion(71,151), 0.0 );
        double ekin=1.e-6;
        TVector3 dir(0,0,1);
        Lu151.setFourMomentumLab(ekin,dir);
        TwoBodyDecay p( Lu151, Ion(1,1), Ion(70,150), std::make_unique<FixedAngularCorrelation>(90,0) );
        p.runProcess();
        auto& d = p.getDaughters();
        TVector3 dir1 = d[0] -> getDirectionLab();
        CHECK_CLOSE( 1., dir1.X(), 1E-3 );
        CHECK_CLOSE( 0., dir1.Y(), 1E-3 );
        CHECK_CLOSE( 0., dir1.Z(), 1E-3 );
    }
    
    TEST_FIXTURE(SetupFixture, OyIsAsExpectedForCompoundDecay) {
        // Process chain
        unique_ptr<ProcessChain> chain = px.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="90" phi="90")
                He4
                Be8
            }
            )");
        // Set beam energy
        Particle& beamPart = chain -> getBeam();
        double ekin = 1.e-6;
        beamPart.setFourMomentumLab(ekin, {0,0,1});
        // Loop over chain and run processes
        for (auto& p : *chain) {
            p.runProcess();
        }
        // Select 8Be->a+a decay
        auto i = begin(*chain);
        ++i;
        // Access alpha particle
        NuclearProcess& proc = *i;
        auto& d = proc.getDaughters();
        auto d1 = d[0];
        TVector3 dir1 = d1 -> getDirectionLab();
        CHECK_CLOSE( 0., dir1.X(), 1E-3 );
        CHECK_CLOSE( 1., dir1.Y(), 1E-3 );
        CHECK_CLOSE( 0., dir1.Z(), 1E-3 );
    }

    TEST_FIXTURE(SetupFixture, Test1BinaryReactionInXZPlane) {
        // Process chain
        unique_ptr<ProcessChain> chain = px.parseString(R"(
            beam: He3
            target: B10
            -> {
                AD: FIXED(theta="45" phi="180")
                p
                C12 Ex: 10840keV
                -> {
                    AD: FIXED(theta="0" phi="0")
                    He4
                    Be8                
                }
            }
            )");

        // Set beam energy
        Particle& beamPart = chain -> getBeam();
        double ekin = 1e-6;
        beamPart.setFourMomentumLab(ekin, {0,0,1});
        // Loop over chain and run processes
        for (auto& p : *chain) {
            p.runProcess();
        }
        // Select 12C->a+8Be decay
        auto i = begin(*chain);
        ++i;
        ++i;
        // Access alpha particle
        NuclearProcess& proc = *i;
        auto& d = proc.getDaughters();
        auto d1 = d[0];
        TVector3 dir1 = d1 -> getDirectionLab();
        CHECK_CLOSE( 1./std::sqrt(2.), dir1.X(), 1E-3 );
        CHECK_CLOSE( 0., dir1.Y(), 1E-3 );
        CHECK_CLOSE( -1./std::sqrt(2.), dir1.Z(), 1E-3 );
    }

    TEST_FIXTURE(SetupFixture, Test2BinaryReactionInXZPlane) {
        // Process chain
        unique_ptr<ProcessChain> chain = px.parseString(R"(
            beam: He3
            target: B10
            -> {
                AD: FIXED(theta="90" phi="180")
                p
                C12 Ex: 10840keV
                -> {
                    AD: FIXED(theta="90" phi="90")
                    He4
                    Be8                
                }
            }
            )");

        // Set beam energy
        Particle& beamPart = chain -> getBeam();
        double ekin = 1e-6;
        beamPart.setFourMomentumLab(ekin, {0,0,1});
        // Loop over chain and run processes
        for (auto& p : *chain) {
            p.runProcess();
        }
        // Select 12C->a+8Be decay
        auto i = begin(*chain);
        ++i;
        ++i;
        // Access alpha particle
        NuclearProcess& proc = *i;
        auto& d = proc.getDaughters();
        auto d1 = d[0];
        TVector3 dir1 = d1 -> getDirectionLab();
        CHECK_CLOSE( 0., dir1.Z(), 1E-3 );
    }

}
