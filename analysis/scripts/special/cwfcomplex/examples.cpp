#include <iostream>
#include <fstream>
#include <string>
#include <complex>

using namespace std;

#define SIGN(a) (((a) < 0) ? (-1) : (1))

const double precision = 1E-10,sqrt_precision = 1E-5;

#include "complex_functions.H"
#include "cwfcomp.cpp"
#include "test_rec_rel.cpp"

int main (void)
{
  string is_it_normalized_str;
  cin>>is_it_normalized_str;

  const bool is_it_normalized = (is_it_normalized_str == "true") ? (true) : (false);
  
  complex<double> l,eta;
  cin>>l>>eta;

  double R;
  cin>>R;

  int Nz,Nl;
  cin>>Nz>>Nl;

  const double step = 2.0/static_cast<double> (Nz);
 
  class Coulomb_wave_functions cwf(is_it_normalized,l,eta),cwf_lp1(is_it_normalized,l+1,eta);

  ofstream out_file("test.output");
  out_file.precision (10);

  if (is_it_normalized) out_file<<"is_it_normalized:true l:"<<l<<" eta:"<<eta<<" R:"<<R<<"  Nz:"<<Nz<<"  Nl:"<<Nl<<endl<<endl;
  if (!is_it_normalized) out_file<<"is_it_normalized:false l:"<<l<<" eta:"<<eta<<" R:"<<R<<" Nz:"<<Nz<<"  Nl:"<<Nl<<endl<<endl;

  for (int i = 0 ; i < Nz ; i++)
  {
    const double x = i*step;
    const complex<double> z = R*exp (complex<double> (0,x*M_PI));

    complex<double> F,dF,Hp,dHp,Hm,dHm,G,dG;
    cwf.F_dF (z,F,dF);
    cwf.G_dG (z,G,dG);
    cwf.H_dH (1,z,Hp,dHp);
    cwf.H_dH (-1,z,Hm,dHm);

    out_file<<"z:"<<z<<endl;
    out_file<<"F:"<<F<<" F':"<<dF<<endl;
    out_file<<"G:"<<G<<" G':"<<dG<<endl;
    out_file<<"H+:"<<Hp<<" H+':"<<dHp<<endl;
    out_file<<"H-:"<<Hm<<" H-':"<<dHm<<endl;
    out_file<<"Wronskian test: "<<Wronskian_test (z,cwf,cwf_lp1)<<endl<<endl;
  }

  complex<double> *const z_tab = new complex<double> [Nz];
  complex<double> * *const F_tab = new complex<double> * [Nz],* *const dF_tab = new complex<double> * [Nz];
  complex<double> * *const G_tab = new complex<double> * [Nz],* *const dG_tab = new complex<double> * [Nz];
  complex<double> * *const Hp_tab = new complex<double> * [Nz],* *const dHp_tab = new complex<double> * [Nz];
  complex<double> * *const Hm_tab = new complex<double> * [Nz],* *const dHm_tab = new complex<double> * [Nz];

  for (int iz = 0 ; iz < Nz ; iz++)
  { 
    z_tab[iz] = R*exp (complex<double> (0,iz*step*M_PI));

    F_tab[iz] = new complex<double> [Nl];
    dF_tab[iz] = new complex<double> [Nl];
    G_tab[iz] = new complex<double> [Nl];
    dG_tab[iz] = new complex<double> [Nl];
    Hp_tab[iz] = new complex<double> [Nl];
    dHp_tab[iz] = new complex<double> [Nl];
    Hm_tab[iz] = new complex<double> [Nl];
    dHm_tab[iz] = new complex<double> [Nl];
  }

  out_file<<endl<<"Recurrence relations results for a table of z values."<<endl;
  cwf_l_tables_recurrence_relations (l,Nl,eta,is_it_normalized,Nz,z_tab,F_tab,dF_tab,G_tab,dG_tab,Hp_tab,dHp_tab,Hm_tab,dHm_tab);

  for (int iz = 0 ; iz < Nz ; iz++)
  {
    out_file<<"z:"<<z_tab[iz]<<endl;
    for (int il = 0 ; il < Nl ; il++) 
    { 
      out_file<<"l[rec]:"<<l+il<<endl;
      out_file<<"F:"<<F_tab[iz][il]<<" F':"<<dF_tab[iz][il]<<endl;
      out_file<<"G:"<<G_tab[iz][il]<<" G':"<<dG_tab[iz][il]<<endl;
      out_file<<"H+:"<<Hp_tab[iz][il]<<" H+':"<<dHp_tab[iz][il]<<endl;
      out_file<<"H-:"<<Hm_tab[iz][il]<<" H-':"<<dHm_tab[iz][il]<<endl<<endl;
    }
  }

  out_file<<endl<<"Recurrence relations results for a single z."<<endl;
  cwf_l_tables_recurrence_relations (l,Nl,eta,is_it_normalized,z_tab[0],F_tab[0],dF_tab[0],G_tab[0],dG_tab[0],Hp_tab[0],dHp_tab[0],Hm_tab[0],dHm_tab[0]);
  out_file<<"z:"<<z_tab[0]<<endl;
  for (int il = 0 ; il < Nl ; il++) 
  { 
    out_file<<"l[rec]:"<<l+il<<endl;
    out_file<<"F:"<<F_tab[0][il]<<" F':"<<dF_tab[0][il]<<endl;
    out_file<<"G:"<<G_tab[0][il]<<" G':"<<dG_tab[0][il]<<endl;
    out_file<<"H+:"<<Hp_tab[0][il]<<" H+':"<<dHp_tab[0][il]<<endl;
    out_file<<"H-:"<<Hm_tab[0][il]<<" H-':"<<dHm_tab[0][il]<<endl<<endl;
  }

  for (int iz = 0 ; iz < Nz ; iz++)
  {
    delete [] F_tab[iz];
    delete [] dF_tab[iz];
    delete [] G_tab[iz];
    delete [] dG_tab[iz];
    delete [] Hp_tab[iz];
    delete [] dHp_tab[iz];
    delete [] Hm_tab[iz];
    delete [] dHm_tab[iz];
  }

  delete [] z_tab;

  delete [] F_tab;
  delete [] dG_tab;
  delete [] G_tab;
  delete [] dF_tab;
  delete [] Hp_tab;
  delete [] dHp_tab;
  delete [] Hm_tab;
  delete [] dHm_tab;

  out_file.close ();
}
