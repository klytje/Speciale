#ifndef KINFITTER_H
#define KINFITTER_H
#include <vector>
#include <functional>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

using namespace std;
using namespace arma;

class KinFitter {
  /**
  * Generic class for kinematic fitting. It needs to be provided with a measurement and
  * its associated covariance matrix. Given any number of constraints it will try to
  * 'fit' the measurement, such that the constraints are met and the deviation between
  * the result and the measurement is minimised. The class relies on the armadillo linear
  * algebra library to do various matrix operations. Before use you should therefore
  * proceed to install armadillo, blas, lapack and other good stuff (which you should
  * want anyway).
  * OBS: This class is not safe to use with the ROOT-interpreter, since proper exception
  *      handling is needed for the matrix invertion algorithms
  */
  private:
    Mat<double> x, y, g;   
    vector<function<double(Mat<double>)>> constraints;  //Constraint functions.
    Mat<double> V;                                      //Covariance matrix.
    Mat<double> G;                                      //Matrix of derivatives.
    double h;        //Small number (used as step for the crude, numerical differentiation).
    int nMax;        //Max no. of iterations allowed for the minimisation.
    double eps;      //Success criterion.
    bool initialised;
    
    void CalculateDerivatives();
    void EvaluateConstraints();

  public:
    KinFitter();
    ~KinFitter();
    
    /**
    * This function is used to add a constraint to the fitter. The fitter will try
    * to find a solution where all the constraints evaluate to zero.
    */
    void AddConstraint(function<double(Mat<double>)>);

    /**
    * When you have added a measurement to the fitter, you are ready to call the Fit()-method,
    * which is the work-horse of the KinFitter. After successful fitting you can retrieve the
    * result with GetResult().
    */
    bool Fit();

    /**
    * Retrieve the result from the kinematic fitting. Should only be called after a measurement
    * has been added with AddMeasurement() and Fit() has been succesfully run.
    */
    vector<double> GetResult();

    /**
    * The fitter uses a crude, homemade, method for differentiation of the multidimensional
    * constraint functions. Basically
    *
    *  (dg/dxi)(x) = (g(...,xi+h/2,...) - g(...,xi-h/2,...))/h
    *
    * If the constraints are well behaved, a small step size should not be a problem.
    */
    void SetDiffStep(double);

    /**
    * Define the success criterion for the minimisation. If the relative improvement
    * in the constructed chi-square is smaller than eps from one iteration to the next,
    * the fitting is regarded as successful and terminates.
    */
    void SetEpsilon(double);

    /**
    * How many iterations you want the minimiser to spend on a single fit.
    */
    void SetMaxIterations(int);

    /**
    * Use this method to give a measurement and its associated covariance-matrix to the
    * fitter. The covariance matrix must be a square matrix with side-length equal to
    * the length of the measurement. Remember to call the Fit()-method before retrieving
    * the result from the fitter.
    */
    void SetMeasurement(vector<double>, Mat<double>);
};
#endif
