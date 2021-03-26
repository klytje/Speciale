
#ifndef coulomb_h
#define coulomb_h

class CoulombFunctions {

    public:
	    CoulombFunctions(double AM1, double AM2, double Z1, double Z2, int Lorb, double R0);
	    virtual ~CoulombFunctions(){}

	    double phi(double e);
	    double penetrability(double e);
	    double shiftFunction(double e);

    private:
	    double am1, am2, z1, z2;
	    int l;
	    double r0;
	    double R, RMAS;

};

#endif 
