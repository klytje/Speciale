script <argument 1> <argument 2> ... <argument n>: description

Dalitz_3D <output/X.root>: creates an interactive 3D Dalitz plot for the given input. This can be used to check if everything looks alright. 

sim_compare <output/X.root> <output/Y.root>: X should be a root file with measured events, while Y should be the simulated events. It creates a Dalitz plot for both inputs and binwise divides each value. Thus the generated figure tells us how the simulation differs from the measurement. 

angular_correlation <j1> <l1>: calculates the angular correlation for 12C --> a + 8Be (ex: 3MeV) --> a + a + a

calibrate_setup <setup/setup.json> <calibration/setup.json>: duplicates the setup file and removes any mention of "calibration". The calibrator tool cannot handle that, while some of the other tools requires it to be present, so this is a hacky solution that keeps everyone happy. 

### ARTICLE FIGURES ###
These scripts have been made specifically for my thesis, and are probably not relevant for anyone else. They should still be functional, however. 

TDC_plot <output/X.root>: compares the TDC values for the four detectors. Note that the locations of the peaks are hardcoded into the script. 
