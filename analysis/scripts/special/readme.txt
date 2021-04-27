script <argument 1> <argument 2> ... <argument n> 
description

Dalitz_3D <output/X.root>
creates an interactive 3D Dalitz plot for the given input. This can be used to check if everything looks alright. 

sim_compare <output/X.root> <output/Y.root>
X should be a root file with measured events, while Y should be the simulated events. It creates a Dalitz plot for both inputs and binwise divides each value. Thus the generated figure tells us how the simulation differs from the measurement. 

angular_correlation <j1> <l1> 
calculates the angular correlation for 12C --> a + 8Be (ex: 3MeV) --> a + a + a

width_converter <reduced width amplitude [keV^(1/2)]> 
calculates the width Gamma corresponding to the reduced width input. Note that there are two implementations: one in python and one in cpp. The python implementation has been tested against a few figures from other papers, while the cpp has not, although it should work just as well. I recommend only using the python implementation anyway, though. 

calibrate_setup <setup/setup.json> <calibration/setup.json> 
duplicates the setup file and removes any mention of "calibration". The calibrator tool cannot handle that, while some of the other tools requires it to be present, so this is a hacky solution that keeps everyone happy. 

ang_compare_single <output path> <nuclear state> <l> <output/X1.root>
plots the theoretical angular correlation function for the given state, and plots it along with a histogram extracted from the data. 

ang_compare <output path> <nuclear state> <l> <4n sim3a files> <4n simX files>
n can be any integer. Essentially the same as ang_compare_single, but it produces figures where it is easy to compare the predictions of the two models. 

theta
plots a Dalitz plot, but with its contents replaced by the angle theta which each point corresponds to. It also produces a plot showing how the angle
depends on the y coordinate. 

interference_merger
needed for sim3a simulations with interference to work, since they do not work out of the box. Made by Morten (like sim3a_i itself).

sim_i_compare
the main script for all sim3a_i simulations. it makes a Dalitz plot and compares some projections of it with the data.

#######################
### ARTICLE FIGURES ###
#######################
These scripts have been made specifically for my thesis, and are probably not relevant for anyone else. They should still be functional, however. 

tdc_plot <output/X1.root> ... <output/XN.root> <output folder> 
compares the TDC values for the four detectors. Note that the locations of the peaks are taken from ../calibrate/plots.cpp (FT & BT), but they can be overwritten locally in any .cpp file.

mul_compare <match/X1.root> ... <match/XN.root> <output folder> 
illustrates how each peak depends on "mul". 

ang_cor_fit <output path> <nuclear state> <l> <output/X.root>
plots a whole bunch of correlation functions along with the histogram extracted from the data. The idea is to reveal which one fits best. I've also
hardcoded a fit for my 0+ data into this file, which should only be produced for 0+ true_events.root data. 

corr_funcs
plots a whole bunch of correlation functions. I had some issues with some of them not being normalized, and this plot helped me realize which.
