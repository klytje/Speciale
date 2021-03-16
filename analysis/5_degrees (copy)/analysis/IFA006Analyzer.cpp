
#include "IFA006Analyzer.h"
#include "TFile.h"
#include "Math/ProbFunc.h"
#include <iostream> // cout
#include <ausa/constants/Mass.h> // cout
#include <ausa/identify/Hit.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/identify/Identification.h>
#include <ausa/util/memory>
#include <TTree.h>

namespace {
    const auto C12_MASS = AUSA::Constants::isotopeMass("C12");
    const auto B11_MASS = AUSA::Constants::isotopeMass("B11");
    const auto BE8_MASS = AUSA::Constants::isotopeMass("Be8");
    const auto PROTON_MASS = AUSA::Constants::isotopeMass("H1");
    const auto ALPHA = AUSA::EnergyLoss::Ion::predefined("He4");

    auto Q = PROTON_MASS + B11_MASS - 3*ALPHA.getMass();
    auto E_SUM = B11_MASS/C12_MASS * 161 + Q;

    auto Z_AXIS = TVector3{0, 0, 1};
}

using namespace std;
using namespace AUSA;
using namespace AUSA::Event;


IFA006Analyzer::IFA006Analyzer()
 : oc()/*, p("TLorentzVector", 3), pCM("TLorentzVector", 3)*/
{
    // make sure everything has been reset
    reset();
}


IFA006Analyzer::~IFA006Analyzer()
{
}


void IFA006Analyzer::setup(std::shared_ptr<Setup> setup)
{
    tree = std::make_unique<TTree>("a", "a");
    tree->Branch("mul", &mul);          //Multiplicitet before reconstruction
	//tree->Branch("exBe8", &exBe8);      //Excitation for Be8 - Todo
    tree->Branch("exC12", &exC12);      //Excitation for C12 - Todo
    tree->Branch("pT", &pTot);          //Total momentum in CM
    tree->Branch("pX", &pX);            //Momentum X
    tree->Branch("pY", &pY);            //Momentum Y
    tree->Branch("pZ", &pZ);            //Momentum Z
    tree->Branch("sV", &sumAng);
    tree->Branch("dV", &devAng);
    tree->Branch("deltaE",&deltaE);

    tree->Branch("zangle", &zangle);

    tree->Branch("N",&N);

    tree->Branch("d",dIndex,"d[3]/I");
	//tree->Branch("FI",FI,"FI[3]/I");
	//tree->Branch("BI",BI,"BI[3]/I");
	//tree->Branch("selected",selected,"selected[3]/I");
	//tree->Branch("time",time,"time[3]/D");
    tree->Branch("eCM",eCM,"eCM[3]/D");
    tree->Branch("eDep",eDep,"eDep[3]/D");
    tree->Branch("eLab",eLab,"eLab[3]/D");
    tree->Branch("thetaCM", thetaCM, "thetaCM[3]/D");
    tree->Branch("thetaLab", thetaLab, "thetaLab[3]/D");
    tree->Branch("phiCM", phiCM, "phiCM[3]/D");
    tree->Branch("phiLab", phiLab, "phiLab[3]/D");
    tree->Branch("vZ", vZ, "vZ[3]/D");
    tree->Branch("vZL", vZL, "vZL[3]/D");
    tree->Branch("exBe8", exBe8, "exBe8[3]/D");
    tree->Branch("prob", prob, "prob[3]/D");



    for (size_t i = 0; i < setup->dssdCount(); ++i) {
        auto name = "loss_" + setup->getDSSD(i)->getName();
        h_loss.push_back(oc.add2D(name, name + ";Energy loss;Energy", 300, 0, 300, 2000, 0, 7E3));
    }
}

void IFA006Analyzer::setScalers(const ScalerOutput &s) {
    scaler = s;
}

void IFA006Analyzer::analyze(const std::vector<PhysicsEvent> &events)
{
    // check that we have at least 1 event interpretation
    if (events.size() == 0) return;


    double eCM[3];
    double vvCM[3];
	double phiCM[3];
    //size_t selected[3];
    TLorentzVector _p[3];

    for (auto& event : events) {
        N = scaler.getOutput("___N___").getValue();

        // we ask for the detection multiplicity of the event
        auto mult = event.getDetectionMultiplicity();

        // increment counter
        if (1 <= mult && mult <= 6) multCounter[mult-1]++;

        if (mult < 2) continue;

        // total energy of beam + target
        const auto pBeam = event.getInitialLorentzVector();
        double ei = pBeam.E();
        auto beamBoost = pBeam.BoostVector();


        TLorentzVector pEjectile;

        mul = 0;
        for (size_t i = 0; i < mult; ++i) {
            auto ion = event.getIon(i);
            if (*ion != ALPHA) continue;
            if (mul==3) break;


            auto& p = event.getLorentzVector(i);

            selected[mul] = i;
            pEjectile += p;
            _p[mul] = p;
            mul++;
        }
        if (!(mul==2 || mul==3)) continue;

        if (mul == 2) {
            _p[2] = pBeam - pEjectile;
            pEjectile+= _p[2];
            selected[2] = -1;
        }

        // total energy of reaction products
        double ef = pEjectile.E();

        // Transform total momentum to CM
        auto pEjectileCM = pEjectile;
        pEjectileCM.Boost(-beamBoost);

        auto p1 = _p[0];
        auto p2 = _p[1];
        this->zangle = p1.Vect().Cross(p2.Vect()).Angle(TVector3(0,0,1))*180/TMath::Pi();


        deltaE = ef-ei;
        exC12 = pEjectileCM.E() - C12_MASS;
        pTot = pEjectileCM.P();
        auto pTag = pTot < 35E3;
        auto eTag = exC12 > 12E3;// && ef-ei<0;

        pX = pEjectile.X();
        pY = pEjectile.Y();
        pZ = pEjectile.Z();
        //


        double max = 0;
        double sum = 0;
        size_t maxI = 0;
        auto primaryTag = false;  // Tag primary alpha in GS decay
        auto beamTag = false;     // Tag if suspected beam particle

        /*
         * Loop over the 3 identified alpha particle
         * Calculates:
         *   CM energy
         *   CM angle
         *
         * This also constructs CM kin curve and fills TTree
         */
        sumAng = 0;
        for (size_t i = 0; i < 3; ++i) {
            //auto n = selected[i];

            //_p[i] = event.getLorentzVector(n);
            auto& p = _p[i];
            //auto id = event.getDetectionEvent()->getIdentification(n);
//            *this->p[i] = p;

            this->vZL[i] = p.Angle(Z_AXIS) * TMath::RadToDeg();

            auto Elab = p.E()-Constants::ALPHA_MASS;
            auto vLab = p.Theta()*180.0/TMath::Pi();
            auto phiLab = p.Phi()*180.0/TMath::Pi();
            p.Boost(-beamBoost);
//            *this->pCM[i] = p;


            if (i > 0)
                sumAng += abs(p.Angle(_p[0].Vect()));
            if (i==2)
                devAng = abs(p.Angle(_p[0].Vect().Cross(_p[1].Vect())));

            auto Ecm = p.E()-Constants::ALPHA_MASS;
            auto vCM = p.Theta()*180.0/TMath::Pi();
            auto phiCM = p.Phi()*180.0/TMath::Pi();
            sum += Ecm;
            eCM[i] = Ecm;
			vvCM[i] = vCM;

            if (Ecm > max) {
                max = Ecm;
                maxI = i;
            }
            if (Ecm > 5600) primaryTag=true;
            //auto eDep = Elab - id->getDeadLayerLoss() - id->getTargetLoss();
            if (Ecm < 360) {
                beamTag = true;
            }

            this->eLab[i] = Elab;
            this->eCM[i] = Ecm;
            //this->eDep[i] = eDep;
            this->thetaCM[i] = vCM;
            this->phiCM[i] = phiCM;
            this->thetaLab[i] = vLab;
            this->phiLab[i] = phiLab;
            //this->dIndex[i] = hit->getDetectorIndex();
            this->vZ[i] = p.Angle(Z_AXIS) * TMath::RadToDeg();
        
			auto pBe8 = pEjectile-p;
        	pBe8.Boost(-pBe8.BoostVector());
        	this->exBe8[i] = pBe8.E() - BE8_MASS;

			this->prob[i] = ROOT::Math::breitwigner_cdf_c(abs(this->exBe8[i]-3030),1513,0)*2;
        }
        sumAng += abs(_p[1].Angle(_p[2].Vect()));

        tree->Fill();
    }

    /*
     * See PhysicsEvent.h for the complete list of get methods !
     */
}


void IFA006Analyzer::terminate()
{
    // write some info
    cout << " Multiplicity counters:" << endl;
    cout << "  * mult-1 : " << multCounter[0] << endl; 
    cout << "  * mult-2 : " << multCounter[1] << endl; 
    cout << "  * mult-3 : " << multCounter[2] << endl; 
    cout << "  * mult-4 : " << multCounter[3] << endl; 
    cout << "  * mult-5 : " << multCounter[4] << endl; 
    cout << "  * mult-6 : " << multCounter[5] << endl;

    gDirectory->WriteTObject(tree.release());
}


void IFA006Analyzer::saveToRootFile(TFileWrapper &file)
{
    oc.write(file);
//    file.write(*tree);
    cout << " Histograms saved to: " << file.get().GetName() << endl;    
}


void IFA006Analyzer::reset()
{
    // reset counter
    for (size_t i=0; i<6; i++) multCounter[i] = 0;
}
