/*********************************************
The R-matrix penetrability and shift function.
Coulomb functions calculated using cwfcomplex
reference: N. Michel, Comp. Phys. Comm. 176 (2007) 232-249
***********************************************/

#include <iostream>

#include "CoulombFunctions.h"
#include "cwfcomp.h"

using namespace std;



CoulombFunctions::CoulombFunctions(double AM1, double AM2, double Z1, double Z2, int Lorb, double R0)
{
  am1 = AM1;
  am2 = AM2;
  z1 = Z1;
  z2 = Z2;
  l = Lorb;
  r0 = R0;
  R = r0*(pow(am1,1./3.) + pow(am2,1./3.));
  RMAS = am1*am2/(am1+am2)*931502.;  
}

double CoulombFunctions::phi(double e)
{
  if (e<=0) return 0;
  double eta, rho, E;
  E = e;
  eta = z1*z2/(137.036)*sqrt(RMAS/(2.*E));
  rho = sqrt(2.*RMAS*E)*R/197329.;
  complex<double> F, dF, G, dG;
  complex<double> lc(l,0.);
  complex<double> zeta(eta,0.);
  complex<double> zrho(rho,0.);
  Coulomb_wave_functions cwf(true,lc,zeta);
  cwf.F_dF (zrho,F,dF);
  cwf.G_dG (zrho,G,dG);
  if (G.real()==0) cout << "Division by zero in CoulombFunctions::phi" << endl;
  return atan((F/G).real());
}  

double CoulombFunctions::penetrability(double e)
{
  if (e<=0) return 0;
  double eta, rho, E;
  E = e;
  eta = z1*z2/(137.036)*sqrt(RMAS/(2.*E));
  rho = sqrt(2.*RMAS*E)*R/197329.;
  complex<double> F, dF, G, dG;
  complex<double> lc(l,0.);
  complex<double> zeta(eta,0.);
  complex<double> zrho(rho,0.);
  Coulomb_wave_functions cwf(true,lc,zeta);
  cwf.F_dF (zrho,F,dF);
  cwf.G_dG (zrho,G,dG);
//  cout << "eta: " << eta << endl;
//  cout << "rho: " << rho << endl;
//  cout << "abs(F*F): " << abs(F*F) << endl;
//  cout << "abs(G*G): " << abs(G*G) << endl;
//  cout << "pen: " << rho/(abs(F*F)+abs(G*G)) << endl;
  return rho/(abs(F*F)+abs(G*G));
}  

double CoulombFunctions::shiftFunction(double e)
{
  //calculate the shiftfunction (in MeV) for an energy e (in MeV) 
  double eta, rho, E;
  E = e;

  complex<double> F, dF, G, dG, Hp, dHp;
  complex<double> lc(l,0.);
  complex<double> i(0.,1.);

  if (E==0)
    {
      /*the Coulomb functions cannot be evaluated for E = 0.
	It is assumed that the functions are continuous at E = 0
	so a mean value is calculated for E = +/- Esmall(keV) 
      */
      double Esmall = 2.;
      eta = z1*z2/(137.036)*sqrt(RMAS/(2.*Esmall));
      rho = sqrt(2.*RMAS*Esmall)*R/197329.;
      complex<double> zeta1(eta,0.);
      complex<double> zeta2(0.,-eta);
      complex<double> zrho1(rho,0.);
      complex<double> zrho2(0.,rho);
      Coulomb_wave_functions cwf1(true,lc,zeta1);
      cwf1.F_dF (zrho1,F,dF);
      cwf1.G_dG (zrho1,G,dG);
      double sh1 = rho*real((F*dF+G*dG)/(abs(F*F)+abs(G*G)));
      Coulomb_wave_functions cwf2(true,lc,zeta2);
      cwf2.H_dH (1.,zrho2,Hp,dHp);
      double sh2 = rho*real(i*dHp/Hp);
      //cout<<sh1<<endl;
      //cout<<sh2<<endl;
      return (sh1+sh2)/2.;
    } 
  if (E>0) 
    {
      eta = z1*z2/(137.036)*sqrt(RMAS/(2.*E));
      rho = sqrt(2.*RMAS*E)*R/197329.;

      complex<double> zeta(eta,0);
      complex<double> zrho(rho,0);
      Coulomb_wave_functions cwf(true,lc,zeta);
      cwf.F_dF (zrho,F,dF);
      cwf.G_dG (zrho,G,dG);
      double sh = rho*real((F*dF+G*dG)/(abs(F*F)+abs(G*G)));
      if (E == 0.2) cout<<sh<<endl;
      return sh;
    }
  else
    {
      /* for negative energy we calculate hankel functions, H+ = G + iF,
	 not Coulomb functions*/
      eta = z1*z2/(137.036)*sqrt(RMAS/(2.*(-E)));
      rho = sqrt(2.*RMAS*(-E))*R/197329.;
      complex<double> zeta(0.,-eta); //imaginary input for E<0
      complex<double> zrho(0.,rho);
      Coulomb_wave_functions cwf(true,lc,zeta);
      cwf.H_dH (1.,zrho,Hp,dHp);
      return rho*real(i*dHp/Hp); //real(complex<double>) returns the real part
    }
}  

