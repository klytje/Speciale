This folder is intended to contain scripts that are not experiment-dependent. Here is a short description of each:
	
	basic: contains c++ scripts for generating figures for single data files. All scripts in this folder are run when "make fig/X" is called
	
	calibrate: contains everything necessary for performing a TDC calibration. This is only meant as a collection of functions, and area thus not executable. Instead you should specify a main method inside your experiment folder. 
		
	special: contains special scripts that are meant to be executed individually, and maybe with weird input. There's another readme in the folder which describes the individual scripts. 
	
	plot_style.cpp is a simple c++ script defining the style of all plots generated.
