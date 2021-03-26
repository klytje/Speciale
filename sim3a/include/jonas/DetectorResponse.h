#ifndef DETECTORRESPONSE_H
#define DETECTORRESPONSE_H
#include <vector>
#include <TF1.h>
#include <TRandom3.h>

using namespace std;

class DetectorResponse {
  /*
  This class emulates the detector response of a 16x16 strip dsssd.
  */
  private:
    vector<TF1> fPhys; //Response functions determining the number of e-h pairs produced in the detector.
    int nBins;         //Length of response function array.
    double Emax;
    TF1 *fFront;       //Response function for F-strips.
    TF1 *fBack;        //Response function for B-strips.
    TF1 *fThres;       //Trigger efficiency (threshold). 100% unless otherwise specified.
    TRandom3 *rGen;    //Random number generator.
    int frontMul, backMul;
    int frontStrips[16], backStrips[16];
    double frontSignals[16], backSignals[16];
    bool frontTriggers[16], backTriggers[16];
    int frontOrder[16], backOrder[16];
    bool stripSort;

    int FindEbin(double);

  public:
    /**
    * The parameters for the constructor are :
    *  -sigPhys: The physical 'resolution' of the detector, straggling of losses through dead layer and e-h pair statistics.
    *  -tau1:    The fast tail of the asymmetric response function.
    *  -tau2:    The slow tail of -------||--------.
    *  -eta:     The weighting between the two tails.
    *  -sigF:    Electronic resolution of the F-strips.
    *  -sigB:    Same for the B-strips.
    *  -Eref:    Reference energy at which the parameters were found.
    *
    * It is known that the peak shapes are energy dependent. As a first guess
    * we can assume the peak asymmetry to be caused solely by the pulse height
    * defect caused by nuclear collisions and use an analytical approximation
    * to the phd to scale the tail parameters as function of energy. If provided
    * with a reference energy (the last argument), the class will modify the peak
    * shapes as function of energy.
    */
    DetectorResponse(double, double ,double, double, double, double, double Eref = -1.);
    ~DetectorResponse();
    
    void SetThreshold(double, double);

    /**
    * This function is called for each event, to see how the detector responds
    * to any number (up to 16) of incoming particles. The arguments are :
    *  -mul: Number of particles.
    *  -fi: Array with information on which front strips were hit.
    *  -bi: Same for back strips.
    *  -energy: The energies available when the particles enter the active volume of the detector.  
    */
    void Generate(int, int *, int *, double *);

    /**
    * 'Getters' for the generated signals.
    */
    int GetFrontMultiplicity();
    int GetBackMultiplicity();
    int GetFrontStrip(int);
    int GetBackStrip(int);
    double GetFrontSignal(int);
    double GetBackSignal(int);
    bool HasFrontTrigger(int);
    bool HasBackTrigger(int);

    /**
    * Should the generated signals be sorted by strip number or not.
    */
    void SortByStrip(bool);
};
#endif
