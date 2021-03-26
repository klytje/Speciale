#ifndef LINESHAPES_H
#define LINESHAPES_H

/**
* This function describes the lineshape of a monoenergetic alpha-line in a Si-detector.
* It is a Gaussian folded with two exponential tails, see for instance 'Collaers and Bortels (1987)'.
* The ordering of the 6 parameters is :
*   -Peak area.
*   -Formal peak centroid (NOT identical to the most probable energy or the mean energy).
*   -Sigma of Gaussian part of the peak.
*   -Fall-off parameter for the dominant (fast) low-energy tail of the peak.
*   -Fall-off parameter for the small (slow) low-energy tail of the peak.
*   -Slow/fast ratio of the tails.
*/
double MonoAlpha(double *, double *);

/**
* This function is a Gaussian peak folded with a single, low-energy exponential tail.
* It can be used for fitting to calibration data, but it also describes implantation
* depth distributions quite well. The ordering of the 4 parameters is :
*   -Peak area.
*   -Formal peak centroid.
*   -Sigma of the Gaussian part of the distribution.
*   -Fall-off parameter for the exponential tail.
*/
double SimpleAlpha(double *, double *);
#endif
