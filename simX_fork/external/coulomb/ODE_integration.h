
#ifndef ODE_integration_h
#define ODE_integration_h

#include "complex_functions.h"



class ODE_integration {

    private:
        std::complex<double> extrapolation_in_zero (const unsigned int n,const std::complex<double> *f) const;
        std::complex<double> F_r_u (const std::complex<double> &r,const std::complex<double> &u) const;
        void integration_Henrici (const unsigned int m,const std::complex<double> &h,
        const std::complex<double> &r0,const std::complex<double> &u0,const std::complex<double> &du0,
        const std::complex<double> &r,std::complex<double> &u,std::complex<double> &du) const;

        const std::complex<double> l,ll_plus_one;  // angular momentum,l(l+1).
        const std::complex<double> two_eta;        // 2.eta, with eta the Sommerfeld parameter.

        unsigned int m_tab[8];                                 // integers used in the extrapolation method.
        double one_over_m_tab[8],interpolation_term_tab[8][8]; // doubles used in the extrapolation method.

    public:
        ODE_integration (const std::complex<double> &l_1,
        const std::complex<double> &two_eta_1);

        void operator() (const std::complex<double> &r0,const std::complex<double> &u0,const std::complex<double> &du0,
        const std::complex<double> &r,std::complex<double> &u,std::complex<double> &du) const;

};




#endif
