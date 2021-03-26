#ifndef complex_functions_h
#define complex_functions_h

#include <math.h>
#include <complex>


// Infinite norm of a complex number.
// ----------------------------------
// It is max(|Re[z]|,|Im[z]|)

inline double inf_norm (const std::complex<double> &z)
{
  return std::max (std::abs (real (z)),std::abs (imag (z)));
}


// Test of finiteness of a complex number
// --------------------------------------
// If real or imaginary parts are finite, true is returned.
// Otherwise, false is returned

inline bool isfinite (const std::complex<double> &z)
{
  const double x = real (z), y = imag (z);

  return (std::isfinite (x) && std::isfinite (y));
}




// Usual operator overloads of complex numbers with integers
// ---------------------------------------------------------
// Recent complex libraries do not accept for example z+n or z==n with n integer, signed or unsigned.
// The operator overload is done here, by simply putting a cast on double to the integer.

inline std::complex<double> operator + (const std::complex<double> &z,const int n)
{
  return (z+static_cast<double> (n));
}

inline std::complex<double> operator - (const std::complex<double> &z,const int n)
{
  return (z-static_cast<double> (n));
}

inline std::complex<double> operator * (const std::complex<double> &z,const int n)
{
  return (z*static_cast<double> (n));
}

inline std::complex<double> operator / (const std::complex<double> &z,const int n)
{
  return (z/static_cast<double> (n));
}

inline std::complex<double> operator + (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline std::complex<double> operator - (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline std::complex<double> operator * (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline std::complex<double> operator / (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)/z);
}








inline std::complex<double> operator + (const std::complex<double> &z,const unsigned int n)
{
  return (z+static_cast<double> (n));
}

inline std::complex<double> operator - (const std::complex<double> &z,const unsigned int n)
{
  return (z-static_cast<double> (n));
}

inline std::complex<double> operator * (const std::complex<double> &z,const unsigned int n)
{
  return (z*static_cast<double> (n));
}

inline std::complex<double> operator / (const std::complex<double> &z,const unsigned int n)
{
  return (z/static_cast<double> (n));
}

inline std::complex<double> operator + (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline std::complex<double> operator - (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline std::complex<double> operator * (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline std::complex<double> operator / (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)/z);
}




inline bool operator == (const std::complex<double> &z,const int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const std::complex<double> &z,const int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) != z);
}







inline bool operator == (const std::complex<double> &z,const unsigned int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const std::complex<double> &z,const unsigned int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) != z);
}



std::complex<double> exp_I_omega_chi_calc (const int omega,const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> sin_chi_calc (const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> log_cut_constant_CFb_calc (const bool is_it_normalized,const int omega,const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> log_cut_constant_CFa_calc (const bool is_it_normalized,const int omega,const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> log_cut_constant_AS_calc (const int omega,const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> log_Cl_eta_calc (const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> sigma_l_calc (const std::complex<double> &l,const std::complex<double> &eta);
std::complex<double> log_Gamma (const std::complex<double> &z);
std::complex<double> log1p (const std::complex<double> &z);
std::complex<double> expm1 (const std::complex<double> &z);


#endif 
