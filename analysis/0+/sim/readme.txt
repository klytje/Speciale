This is the directory for simulations. 

configs: contains the basic descriptions of the simulations we want to perform. The extension denotes what simulator will be used to perform it.
	.simX: The standard simX library will be used. 
	.balamuth: The standard sim3a based on the balamuth equation will be used.
	.balamuth_i: (UNSTABLE!) The sim3a modified by a previous student to support interference will be used. 

Once a configuration has been saved into the configs folder, "make fig/X" will perform all of the necessary operations and output a set of figures in ../figures/. Note that makefile is not smart enough to figure out which simulator to use if multiple configurations have the same name (eg. 0+.simX, 0+.balamuth, 0+.balamuth_i), so a suffix is required. I've decided to use no prefix for sim3a, _simX for simX, and _i for sim3a_i, which are used in the makefile logic.

data: contains the simulated data from simX.

raw: contains the raw output from sim3a, which will be used in simX. 

temp: contains some temporary files for my thesis-specific simulations. Everything in it will be generated on the fly, and can thus safely be deleted. 
