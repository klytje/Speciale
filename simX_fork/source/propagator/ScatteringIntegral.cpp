
#include "simX/propagator/ScatteringIntegral.h"
#include <ausa/util/memory>
#include <cmath>

using namespace std;
using namespace simX;
using namespace simX::propagator;


// z range:     +0.002 to +40
// ln(z) range: -6.215 to +3.689 


ScatteringIntegral::ScatteringIntegral(double zmin, double zmax)
 : zmin(zmin),
   zmax(zmax),
   delta(1e-12)   
{
    // load tabulation of f(z)
    loadTabulation();

    // interpolate f(ln(z))    
    fzInter = make_unique<ROOT::Math::Interpolator> (ROOT::Math::Interpolation::kCSPLINE);
  	Double_t xlh[15];
	for ( Int_t i=0; i<15; ++i ) {
		xlh[i] = std::log(zlh[i]);
	}
	fzInter -> SetData(15, xlh, flh);

	// Define J(z) as in Eq. (5)
	double xmin = std::log(1.e-5);
	double xmax = std::log(1.e5);
	TF1 * dJdx = new TF1("integrand", this, &ScatteringIntegral::finteq5, xmin, xmax, 0);

	// Tabulate J(z)
	int ntab = 1E3;
	double z0 = 0.05;
	double x0 = std::log(z0);
	double ztab[ntab], Jztab[ntab];
	xmin = std::log(zmin);
	xmax = std::log(zmax);
	for (int i=0; i<ntab; ++i) 
	{   
		double x = xmin + (float)i / (float)(ntab-1) * (xmax-xmin);
		ztab[i]  = std::exp(x);
		Jztab[i] = dJdx -> Integral(x0, x);
	}

	// Interpolate J(z)
    jzInter = make_unique<ROOT::Math::Interpolator> (ROOT::Math::Interpolation::kCSPLINE);
	jzInter -> SetData(ntab, ztab, Jztab);

	// Define TF1 for J(z)
	Jz = new TF1("Jz", this, &ScatteringIntegral::fjz, zmin, zmax, 0);
}


ScatteringIntegral::~ScatteringIntegral()
{
    // do nothing
}


double ScatteringIntegral::eval(double z)
{
    return Jz -> Eval(z);
}


double ScatteringIntegral::invert(double j, double zlow, double zhigh)
{
    double zl = max(zlow, zmin);
    double zh = min(zhigh, zmax);
    
    return Jz -> GetX(j, zl, zh);
}


void ScatteringIntegral::loadTabulation() 
{
    zlh[0] = 0.002;
    zlh[1] = 0.004;
    zlh[2] = 0.01;
    zlh[3] = 0.02;
    zlh[4] = 0.04;
    zlh[5] = 0.10;
    zlh[6] = 0.15;
    zlh[7] = 0.20;
    zlh[8] = 0.40;
    zlh[9] = 1;
    zlh[10] = 2;
    zlh[11] = 4;
    zlh[12] = 10;
    zlh[13] = 20;
    zlh[14] = 40;

    flh[0] = 0.162;
    flh[1] = 0.209;
    flh[2] = 0.280;
    flh[3] = 0.334;
    flh[4] = 0.383;
    flh[5] = 0.431;
    flh[6] = 0.435;
    flh[7] = 0.428;
    flh[8] = 0.385;
    flh[9] = 0.275;
    flh[10] = 0.184;
    flh[11] = 0.107;
    flh[12] = 0.050;
    flh[13] = 0.025;
    flh[14] = 0.0125;
}


Double_t ScatteringIntegral::finteq5(Double_t *x, Double_t *par ) 
{
    Double_t f = 0;
    Double_t z = std::exp(x[0]);

    if (z > zlh[0] && z < zlh[14])  f = fzInter -> Eval(x[0]);
    else if (z <= zlh[0])           f = 0;       // bad approximation, should be improved
    else if (z >= zlh[14])          f = 0.5/z;   // pure rutherford

    return f/z;
}


Double_t ScatteringIntegral::fjz(Double_t *z, Double_t *par) 
{
    Double_t jz = 0;

    if (z[0] > zmin && z[0] < zmax) jz = jzInter -> Eval(z[0]);
    else if (z[0] <= zmin)           jz = jzInter -> Eval(zmin + delta);
    else if (z[0] >= zmax)           jz = jzInter -> Eval(zmax - delta);
 
    return jz;
}

