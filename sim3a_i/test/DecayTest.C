#include <TRandom3.h>
#include <TH2D.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <iostream>
#include <ausa/constants/Mass.h>
#include "./include/Decay.h"

using namespace std;
using namespace AUSA;
using namespace Constants;

TH2D * test()
{
  Decay *event = new Decay();
  event->SetType(2);
  event->DoRecoil(false);
  double energy = 3*isotopeMass(2,4) + 2000. - isotopeMass(6,12);
  event->SetEnergy(energy);
 
  double alphaMass = isotopeMass(2,4);
  TH2D *test = new TH2D("test","",140,0,1.4,140,-0.7,0.7);
  TRandom3 *rGen = new TRandom3(0);

  TH1I *h = new TH1I("h","",1000,0,1);

  for(int i=0; i<1e2; i++){
    double weight = event->Generate();
    double max = event->GetMaxWeight();
    double num = rGen->Uniform(0,max);

    h->Fill(weight);
    if(weight < num) continue;

    double EA1 = event->GetCmEnergy(0);
    double EA2 = event->GetCmEnergy(1);
    double EA3 = event->GetCmEnergy(2);

    double ea1 = event->GetProduct(0)->Energy()*1e6 - alphaMass;
    double ea2 = event->GetProduct(1)->Energy()*1e6 - alphaMass;
    double ea3 = event->GetProduct(2)->Energy()*1e6 - alphaMass;

    TLorentzVector total = *(event->GetProduct(0)) + *(event->GetProduct(1)) + *(event->GetProduct(2));
    double Etot = total.Energy()*1e6 - 3*alphaMass;

    /*
    cout << "E1 = " << EA1 << endl;
    cout << "E2 = " << EA2 << endl;
    cout << "E3 = " << EA3 << endl;
    cout << "Etot = " << EA1 + EA2 + EA3 << endl;
    cout << "In lab system: " << endl;
    cout << "E1 = " << ea1 << endl;
    cout << "E2 = " << ea2 << endl;
    cout << "E3 = " << ea3 << endl;
    cout << "Etot = " << Etot << endl;
    cout << endl;
    */

    //cout << weight << endl;

    int alphaOrder[3];
    double alphaEnergies[3] = {EA1, EA2, EA3};
    TMath::Sort(3,alphaEnergies,alphaOrder,true);
    double E1 = alphaEnergies[alphaOrder[0]];
    double E2 = alphaEnergies[alphaOrder[1]];
    double E3 = alphaEnergies[alphaOrder[2]];
    double totalEnergy = E1 + E2 + E3;
    double x = Sqrt(3) * (E1 - E3) / totalEnergy;
    double y = (2*E2 - E1 - E3) / totalEnergy;
    test->Fill(x,y);
  }

  h->Draw();

  return test;
}
