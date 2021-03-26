#include <iostream>
#include "./include/Nucleus.h"
#include "./include/Level.h"
#include "./include/ParticlePair.h"
#include "./include/SeqDecay.h"

using namespace std;

void test(int N = 1e5)
{
  double Q1 = (12710. - (7366.6 - 91.8)); //Available energy (keV) for 3alpha-breakup.

  Nucleus *C12 = new Nucleus(12,6,0,1);
  Nucleus *Be8 = new Nucleus(8,4,0,1);
  Nucleus *He4 = new Nucleus(4,2,0,1);

  Level *mother = new Level(12710.,1,1);
  ParticlePair *pair1 = new ParticlePair(Be8,He4);
  int pChannelID = mother->AddParticleChannel(pair1,2,30.,6.7);
  int pLevelID = C12->AddLevel(mother);

  Level *intermediate = new Level(3037.,2,1);
  ParticlePair *pair2 = new ParticlePair(He4,He4);
  int sChannelID = intermediate->AddParticleChannel(pair2,2,32.787,4.5);
  int sLevelID = Be8->AddLevel(intermediate);

  SeqDecay *decay = new SeqDecay();
  decay->SetPrimary(C12,pLevelID,pChannelID);
  decay->SetSecondary(Be8,sLevelID,sChannelID);
  double Q2 = decay->Generate();

  cout << "Q1 = " << Q1 << ",  Q2 = " << Q2 << endl;

  /*
  double alphaMass = 4He->M() / 1e6;
  double beMass = 8Be->M() / 1.e6;
  double masses[3] = { alphaMass, alphaMass, alphaMass};
  TLorentzVector mother(0.0,0.0,0.0,3*alphaMass + Q/1.e6);

  //Decay from 12C(1+) state through 8Be(2+).
//  int l1 = 2; //First emitted alpha. 
//  int l2 = 2; //Second...

//  int ja = 1; //Spin of 12C.
//  int jb = 2; //Spin of 8Be

//  double r1 = 6.7;    //Channel radius of 12C-decay in fm.
//  double r2 = 4.5;    //Channel radius of 8Be-decay in fm.
//  double Eb = 3129.;  //Excited level in 8Be (measured from 2alpha-threshold).
//  double g2a = 1000.; //Arbitrary number, only penetrability-independence important.
//  double g2b = 1075.; //Squared reduced width in keV

//  ParticleChannel c1(alphaMass*1e6,beMass*1e6,2,4,l1,r1);
//  ParticleChannel c2(alphaMass*1e6,alphaMass*1e6,2,2,l2,r2);
//  ClebschGordan cg(l1,jb);
  
  TH2I *h = new TH2I("h","",140,0,1.4,140,-0.7,0.7);
  TH1I *g = new TH1I("g","",2000,0,2);
  TH1I *f = new TH1I("f","",1200,0,6000);
  TRandom3 rGen(0);

  for(int k=0; k<N; k++){
    TGenPhaseSpace event;
    event.SetDecay(mother, 3, masses);
    double weight = event.Generate();
    double alphaEnergies[3];
    alphaEnergies[0] = (event.GetDecay(0)->Energy() - alphaMass) * 1e6;
    alphaEnergies[1] = (event.GetDecay(1)->Energy() - alphaMass) * 1e6;
    alphaEnergies[2] = (event.GetDecay(2)->Energy() - alphaMass) * 1e6;
    int alphaOrder[3];
    Sort(3,alphaEnergies,alphaOrder,true);
    double  e1 = alphaEnergies[alphaOrder[0]];
    double  e2 = alphaEnergies[alphaOrder[1]];
    double  e3 = alphaEnergies[alphaOrder[2]];
    double x = Sqrt(3) * (e1 - e3) / (e1 + e2 + e3);
    double y = (2*e2 - e1 - e3) / (e1 + e2 + e3); 
    
    double prob = 0;
    //Initial spin direction {-1,0,1}, to be averaged over.
    for(int ma = -ja; ma<=ja ;ma++){
      complex<double> amplitude(0,0);
      for(int j=0; j<3; j++){
        TLorentzVector *a1 = event.GetDecay(j);
        TLorentzVector *a2 = event.GetDecay((j+1)%3);
        TLorentzVector *a3 = event.GetDecay((j+2)%3);
    
        TLorentzVector alpha1 = *a1;
        TLorentzVector alpha2 = *a2;
        TLorentzVector alpha3 = *a3;

        double E1 = (alpha1.Energy() - alphaMass) * 1e6;
        TLorentzVector rcm = alpha2 + alpha3;     //Recoil center of mass system.

        double theta1 = alpha1.Theta();
        double phi1 = alpha1.Phi();
        alpha2.Boost(-rcm.BoostVector());
        double theta2 = alpha2.Theta();
        double phi2 = alpha2.Phi();

        rcm.Boost(-rcm.BoostVector());
        double E23 = (rcm.Energy() - 2 * alphaMass) * 1e6;  //Relative energy of alpha2 and alpha3.
        double E = Q - E23;                                 //Relative energy of alpha1 and recoil.

        double Gamma1 = 2 * c1.Penetrability(E) * g2a;
        double Gamma2 = 2 * c2.Penetrability(E23) * g2b;
        double phase1 = c1.CoulombShift(E) - c1.HardSphereShift(E);
        double phase2 = c2.CoulombShift(E23) - c2.HardSphereShift(E23);
  
        //We start calculating...
        complex<double> i(0,1);
        complex<double> a(0,0);
        a += sqrt(Gamma1 * Gamma2 / sqrt(E1 * E23)) * exp(i * phase1) * exp(i * phase2);
        a /= Eb - g2b * (c2.ShiftFunction(E23) - c2.ShiftFunction(Eb)) - E23 - 0.5 * i * Gamma2;
        complex<double> sum(0,0);
        for(int mb = -jb ; mb <= jb; mb ++){
          complex<double> Cmmj(cg.Coefficient(ma-mb,mb,ja),0);
          complex<double> Ylm1(spherical_harmonic(l1,ma-mb,theta1,phi1));
          complex<double> Ylm2(spherical_harmonic(l2,mb,theta2,phi2));
          sum += Cmmj * Ylm1 * Ylm2;
        }
        amplitude += a * sum;
      }
      prob += norm(amplitude);
    }
    prob *= 11440.;   //Apparent maximum probability;
    g->Fill(prob);
    weight /= 0.5;  //Maximum weight.
    double totalWeight = weight * prob;
    double num = rGen.Uniform(0,1);
    if(num < totalWeight){
      h->Fill(x,y);
      f->Fill(e1); f->Fill(e2); f->Fill(e3);
    }
  }
  g->Draw();
  return h;
  */
}
