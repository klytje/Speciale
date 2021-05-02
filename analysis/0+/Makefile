# AUSAlib executables. If globally installed, replace with "./"
sorter := ~/tools/sorter/build/Sorter
calibrator := ~/tools/calibrator/build/Calibrator
identifier := ~/tools/identifier/build/Identifier

.DEFAULT_GOAL := help
.SECONDARY:
.PHONY: clean clean/build clean/sim_widths

help: 
	@ echo "Main targets: (consider using the -jX flag if data/ contains many files)"
	@ echo "\"build\": Prepares all executables for compilation. Must be run once to initialize the project. \n"
	@ echo "\"fig/X\": Generates the figures for data source X. data/X.root, output/X.root, sim/configs/X.bahamut, or sim/configs/X.simX must exist. \n"
	@ echo "\"tdc_calibration\": Takes all files in data/ and performs a TDC calibration on them. Outputs output/true_events.root, output/reconstructed_events.root, and output/all_events.root which can be used as targets for another call to fig/X. \n"
	@ echo "\"clean/X\": Removes all generated data and figures associated with data source X, leaving only the raw source in data/ or sim/config/ . \n"
	@ echo "\"clean\": Removes all generated files, leaving only the raw data sources. \n"
	@ echo "There are also further figures available; take a look at the Makefile section \"figures\" for an overview. \n"
	@ echo "Note that for any AUSAlib tool to work, you must update their executable paths in the first few lines of the Makefile. \n"

# removes everything, leaving only the raw data files in data/
clean: clean/build clean/sim_widths
	rm -f sim/raw/*
	rm -f sim/data_i/*
	rm -f sim/data/*
	rm -f match/*
	rm -f merged/*
	rm -f id/*
	rm -f output/*
	rm -f -r figures/*
	rm -f -r analysis/figures/*

# removes all files associated with %, except for the raw data in data/
# note: also removes any simulated data named %.root
clean/%: 
	@ rm -f sim/raw/$(*F).dat
	@ rm -f sim/data_i/$(*F).root
	@ rm -f sim/data/$(*F).root
	@ rm -f match/$(*F)m.root
	@ rm -f id/$(*F)mi.root
	@ rm -f output/$(*F).root
	@ rm -f merged/$(*F).root
	@ rm -f -r figures/$(*F)
#	@ rm -f data/$(*F).root

#############
### CMAKE ###
#############
../../simX_fork/build/Makefile: ../../simX_fork/CMakeLists.txt
	mkdir -p ../../simX_fork/build
	cd ../../simX_fork/build; cmake ../
	
../../sim3a/build/Makefile: ../../sim3a/CMakeLists.txt
	mkdir -p ../../sim3a/build
	cd ../../sim3a/build; cmake ../

../../sim3a_i/build/Makefile: ../../sim3a_i/CMakeLists.txt
	mkdir -p ../../sim3a_i/build
	cd ../../sim3a_i/build; cmake ../

../scripts/basic/build/Makefile: ../scripts/basic/CMakeLists.txt
	mkdir -p ../scripts/basic/build
	cd ../scripts/basic/build; cmake ../

../scripts/special/build/Makefile: ../scripts/special/CMakeLists.txt
	mkdir -p ../scripts/special/build
	cd ../scripts/special/build; cmake ../

analysis/build/Makefile: analysis/CMakeLists.txt
	mkdir -p analysis/build
	cd analysis/build; cmake ../

# this variable keeps track of changes in all CMakeLists.txt and rebuilds if necessary
# every time a CMakeLists is modified, run this target
build: ../../simX_fork/build/Makefile ../../sim3a/build/Makefile ../../sim3a_i/build/Makefile ../scripts/basic/build/Makefile ../scripts/special/build/Makefile analysis/build/Makefile

####################
### DEPENDENCIES ###
####################
simX := ../../simX_fork/build/simX
balamuth := ../../sim3a/build/balamuth
balamuth_i := ../../sim3a_i/build/balamuth_i

# if any of the dependent executables doesn't exist, we compile them
../../simX_fork/build/simX_fork/build/simX:
	make -C ../../simX_fork/build/
	
../../sim3a/build/balamuth: ../../sim3a/balamuth/*
	make -C ../../sim3a/build/

../../sim3a_i/build/balamuth: ../../sim3a_i/balamuth/*
	make -C ../../sim3a_i/build/

# removes all executables
clean/build: 
	@ rm -f -r ../figure_scripts/build/*
	@ rm -f -r ../special_figure_scripts/build/*
	@ rm -f -r ../calibrate/build/*
	@ rm -f -r analysis/build/*

#######################
### tdc_calibration ###
#######################
calibrate_rec := analysis/build/tdc_calibration_reconstructed
calibrate_true := analysis/build/tdc_calibration_true

# this target generates a whole bunch of figures in analysis/figures/, which you should check for errors. As long as none of the plots are empty, it should be ok.
tdc_calibration: $(calibrate_rec) $(calibrate_true) $(addprefix output/,$(shell ls data))
	$(calibrate_true) $(addprefix output/,$(shell ls data))	
	$(calibrate_rec) $(addprefix output/,$(shell ls data))
	@ echo "Creating output/all_events.root ..."
	@ hadd -f output/all_events.root output/true_events.root output/reconstructed_events.root

### compile the calibration programs if necessary ###
analysis/build/tdc_calibration_reconstructed: $(build) analysis/tdc_calibration_reconstructed.cpp ../scripts/calibrate/*
	make -C analysis/build/ tdc_calibration_reconstructed

analysis/build/tdc_calibration_true: $(build) analysis/tdc_calibration_true.cpp ../scripts/calibrate/*
	make -C analysis/build/ tdc_calibration_true
	
#####################
### DEBUG TARGETS ###
#####################
debug_match: $(calibrate_rec)
	../calibrate/build/mul_compare $(addprefix match/,$(shell ls match))

test: $(calibrate_rec) $(calibrate_true)
	$(calibrate_rec) output/2455_000.root output/2455_001.root output/2455_002.root

test_all: $(calibrate_rec) $(calibrate_true) $(addprefix output/,$(shell ls data))
	$(calibrate_rec) $(addprefix output/,$(shell ls data))

temp: calibration/setup.json
	$(calibrator) -s $< -k setup/source.json -d SU -o calibration/ calibration/SU/2478_000.root

#######################
### REAL ROOT FILES ###
#######################
### calibration ###
# performs a calibration of detector %, using all root files present in the folder calibration/%/
# it uses a small script to duplicate the setup file from setup/setup.json, while removing any mention of "calibration: ", 
# since this breaks the calibrator, but is required for some of the other tools. 
setup/%.cal: ../scripts/special/calibrate_setup.py setup/setup.json
	python3 $< setup/setup.json calibration/setup.json
	$(calibrator) -s calibration/setup.json -k setup/source.json -d $* -o calibration/ --eloss --gui calibration/$*/*.root 
	@ cp calibration/*.cal setup/
	@ rm calibration/setup.json

# simple target to calibrate all detectors
calibrate: setup/SU.cal setup/SD.cal setup/Det1.cal setup/Det2.cal;

### matching ###
match/%m.root : data/%.root matcher.json setup/setup.json calibrate
	$(sorter) $< -m matcher.json

# this target allows us to use simulated data from this point on
match/%m.root : sim/data/%.root matcher.json setup/setup.json calibrate
	$(sorter) $< -m matcher.json

### identification ###
collectors = --collectors SU:EA:TB,SD:EA:TB,Det1:EA:TF,Det2:EA:TF
target = -t setup/target.json
ions = -i "a"
margin = -m 500keV
id/%i.root : match/%.root
	$(identifier) $(target) $(ions) $(margin) $(collectors) $^ -o id/$*i.root

### analysis ###
event_analysis := analysis/build/EventAnalysis

analysis/build/EventAnalysis: analysis/StandardAnalyzer.cpp analysis/main.cpp
	make -C analysis/build/ EventAnalysis

# if the output is based on a simulation with sim3a_i, we need to perform a special merge operation on it as well
output/%.root : id/%mi.root sim/configs/%.balamuth_i $(event_analysis)
	$(event_analysis) output/$*.root id/$*mi.root match/$*m.root
	echo sim/raw/$*.root output/$*.root | root -l ../scripts/special/interference_merger.C

# otherwise we just run the standard analyzer on it and do nothing else
output/%.root : match/%m.root id/%mi.root $(event_analysis)
	$(event_analysis) output/$*.root id/$*mi.root match/$*m.root

###############
### FIGURES ###
###############
# the following line allows us to use automatic variables for prerequisites from this point on
.SECONDEXPANSION:

# from here on I'll make heavy use of a few somewhat exotic automatic variables.
# example: pattern: figures/%.pdf, matched by figures/abc/xyz.pdf
# @D --> figures/abc
# *D --> abc
# @F --> xyz.pdf
# *F --> xyz
# $* --> abc/xyz

# compiles a given figure program
../scripts/basic/build/%: ../scripts/basic/%.cpp ../scripts/calibrate/*
	make -C ../scripts/basic/build/ $*

# generates a single pdf file for the data source %
figures/%.pdf: output/$$(*D).root ../scripts/basic/build/$$(*F)
	@ echo "Building "$@"..."
	@ mkdir -p figures/$(*D)
	../scripts/basic/build/$(*F) $< $(shell pwd)/figures/$(*D)/

# contains a list of all possible pdf files
figure_pdf = $(addsuffix .pdf,$(basename $(shell find ../scripts/basic/ -maxdepth 1 -name "*.cpp" -printf "%P ")))

# concatenates all raw data in data/ and makes all combined figures for them
figs: $(addprefix figures/data/,$(figure_pdf))

# fig/X generates all figures for X. Either data/X.root, ../sim/X.simX, ../sim/X.bahamut must exist
# this command will simulate and/or perform all necessary analysis of the data, and is thus a nice
# point of entry to this whole project 
fig/% : $(addprefix figures/%/,$(figure_pdf)) output/%.root
	@ echo "Built figures/"$*

### special figures ###
# compiles a given "special" figure program
../scripts/special/build/%: ../scripts/special/%.cpp ../scripts/utility.cpp
	make -C ../scripts/special/build/ $*

# generates the figures for a sim3a_i simulation (due to the way they are generated, the standard figures does not work for them)
fig/%_i: 
	@ mkdir -p figures/$*_i
	$(MAKE) sim_compare/$*_i

# creates an interactive 3D Dalitz plot of %. mainly used for debug
Dalitz_3D/%: output/%.root ../scripts/special/build/Dalitz_3D
	@ echo "Creating interactive 3D Dalitz plot..."
	../scripts/special/build/Dalitz_3D $<

# compare a simulation with the data. Generates a variety of projection plot for easy comparison
sim_compare/%: ../scripts/special/build/sim_compare output/%.root output/true_events.root
	$< $(shell pwd)/figures/$*/ output/$*.root output/true_events.root

# plots before / after of the TDC values
tdc_plots: ../scripts/special/build/tdc_plot $(addprefix output/,$(shell ls data))
	@ echo "Creating figures/article/tdc_plot_before.pdf and figures/article/tdc_plot_after.pdf"
	$< $(addprefix output/,$(shell ls data)) $(shell pwd)/figures/article/

# illustrates how the peaks depends on the multiplicity
depend_mul: ../scripts/special/build/depend_mul $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data))))
	@ echo "Creating figures/article/depend_mul.pdf"
	$< $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data)))) $(shell pwd)/figures/article/

# illustrates how the peaks depends on the energy
depend_E: ../scripts/special/build/depend_E $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data))))
	@ echo "Creating figures/article/depend_E.pdf"
	$< $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data)))) $(shell pwd)/figures/article/
#	$< match/2455_000m.root match/2455_001m.root match/2455_002m.root $(shell pwd)/figures/article/

# shows how theta' changes throughout the dalitz plot
theta: ../scripts/special/build/theta
	$< $(shell pwd)/figures/article/

# calculates theta based on the data
theta/%: ../scripts/special/build/theta_data output/%.root
	$< $(shell pwd)/figures/$*/ output/$*.root

# shows a bunch of different correlation functions. This was mainly to verify that they were correctly normalized
corr_funcs: ../scripts/special/build/corr_funcs
	$< $(shell pwd)/figures/article/

# note that the next few targets depends on sim/temp/converted_widths.txt. Due to limitations of make (and/or my ability), this file must be created manually before any target below is executed. If you mess it up, simply delete the sim/temp folder, or make the target clean/sim_widths.
# extra information for the next couple of targets
state := 0+
l := 2

# compare multiple Dalitz plots from simX and sim3a with equivalent widths
depend_width: ../scripts/special/build/width_compare sim/temp/output
	$< $(shell pwd)/figures/article/ $(addprefix output/0+_, $(addsuffix .root, $(shell awk '{print $$1}' sim/temp/converted_widths.txt | awk -vORS=" " 'NR>1 {print}'))) $(addprefix output/0+_, $(addsuffix .root, $(shell awk '{print $$2}' sim/temp/converted_widths.txt | awk -vORS=" " 'NR>1 {print}')))

# compare the angular correlation implicit in sim3a with the theoretical prediction for a single data source. this is useful to check if everything looks alright before doing the similar ang_compare for all data sources. 
ang_compare/%: ../scripts/special/build/ang_compare_single output/%.root
	$< $(shell pwd)/figures/$*/ $(state) $(l) output/$*.root

# same as ang_compare, but plots a whole bunch of correlation functions to determine which one fits best
# should really only be used with actual data; simulated data should match ang_compare very well (otherwise something is really wrong)
ang_cor_fit/%: ../scripts/special/build/ang_cor_fit output/%.root
	$< $(shell pwd)/figures/$*/ $(state) $(l) output/$*.root

# same as ang_compare_single, but for all files simulated with width_compare
ang_compare: ../scripts/special/build/ang_compare sim/temp/output
	$< $(shell pwd)/figures/article/ $(state) $(l) $(addprefix output/0+_, $(addsuffix .root, $(shell awk '{print $$1}' sim/temp/converted_widths.txt | awk -vORS=" " 'NR>1 {print}'))) $(addprefix output/0+_, $(addsuffix .root, $(shell awk '{print $$2}' sim/temp/converted_widths.txt | awk -vORS=" " 'NR>1 {print}')))

###################
### SIMULATIONS ###
###################
# number of events generated by all simulations
N := 1000000

# generate the simX configuration file required for simX, with data from sim3a
# this is generated on the fly whenever it is needed, and is deleted afterwards
sim/%.simX: sim/raw/%.dat
	@ echo "beam: p" > $@
	@ echo "target: B11" >> $@
	@ echo "-> {" >> $@
	@ echo "    generator: PREGEN(file=\""$(shell pwd)/$<"\" unit=\"keV\")" >> $@
	@ echo "    a" >> $@
	@ echo "    a" >> $@
	@ echo "    a" >> $@
	@ echo "}" >> $@

# generate the simX configuration file required for simX, with data from sim3a with interference
sim/%.simX: sim/raw/%.root
	@ echo "beam: p" > $@
	@ echo "target: B11" >> $@
	@ echo "-> {" >> $@
	@ echo "    generator: PREGEN(file=\""$(shell pwd)/$<"\" unit=\"keV\")" >> $@
	@ echo "    a" >> $@
	@ echo "    a" >> $@
	@ echo "    a" >> $@
	@ echo "}" >> $@

# generate the raw event data from sim3a
sim/raw/%.dat: sim/configs/%.balamuth $(balamuth)
	$(balamuth) $(shell cat $<) -N $(N) -o $(shell pwd)/sim/raw/$*.dat	

# simulate using balamuth with interference
sim/raw/%.root: sim/configs/%.balamuth_i $(balamuth_i)
	$(balamuth_i) $(shell cat $<) -N $(N) -o $(shell pwd)/$@

# simulate the detector setup with simX, using the generated events from sim3a (possibly with interference)
sim/data/%.root: sim/%.simX $(simX)
	$(simX) -c sim/config.json -r sim/$*.simX -o sim/data/$*.root -N $(N)
	@ rm $<

# simulate everything solely with simX
sim/data/%.root: sim/configs/%.simX $(simX)
	$(simX) -c sim/config.json -r $< -o sim/data/$*.root -N $(N)

###################################
### thesis specific simulations ###
###################################
# I sometimes create new files (i.e. sim/temp/setup) which are not used for anything in this section. The idea is that they are just shortcuts for longer operations, and Makefile needs to know when such a target was last executed. Creating a new file is the perfect way to do just that.

# calculates a few different values of the width Gamma based on the reduced width amplitude gamma
sim/temp/converted_widths.txt: ../scripts/special/width_converter.py
	@ mkdir -p sim/temp
	@ echo "gamma Gamma" > $@
	@ for g in 1 2 3 4 5 10 15 20 25 30 35 40; do python3 $< $$g -l >> $@; done

sim/temp/0+_%.balamuth:
	@ echo "--JC 0+ --JBe 2+ -L 2 --ExC 17.76MeV --ExBe 3030keV --gBe $*keV" > $@

# generate the raw event data from sim3a
sim/raw/%.dat: sim/temp/%.balamuth $(balamuth)
	$(balamuth) $(shell cat $<) -N $(N) -o $(shell pwd)/$@

sim/temp/0+_%.simX:
	@ echo "beam: p" > $@
	@ echo "target: B11" >> $@
	@ echo "-> {" >> $@
	@ echo "    a" >> $@
	@ echo "    Be8 Ex: 3.03MeV G0: $*MeV" >> $@
	@ echo "    -> {" >> $@ 
	@ echo "        AD: TF1(phi=\"1\" theta=\"sin(x)*(11.25*cos(x)^4 - 7.5*cos(x)^2 + 1.25)\")" >> $@
	@ echo "        a" >> $@
	@ echo "        a" >> $@
	@ echo "       }" >> $@
	@ echo "   }" >> $@
	@ echo "}" >> $@
	
# this is just a bit of logic that creates all the simX & sim3a configuration files in the sim/temp folder. 
setup: sim/temp/converted_widths.txt $(addprefix sim/temp/0+_, $(addsuffix .simX, $(shell awk '{print $$2}' sim/temp/converted_widths.txt | awk -vORS=" " 'NR>1 {print}'))) $(addprefix sim/temp/0+_, $(addsuffix .balamuth, $(shell awk '{print $$1}' sim/temp/converted_widths.txt | awk -vORS=" " 'NR>1 {print}')));

# simulate everything solely with simX
sim/data/%.root : sim/temp/%.simX $(simX)
	$(simX) -c sim/config.json -r $< -o $@ -N $(N)

# same explanation as above, except this time it prepends "output/0+_" and appends ".root", such that the entire AUSAlib pipeline is executed on all simulations. note that this takes *a lot* of time to do - I highly recommend using the -jX flag for this!
sim/temp/output: setup $(addprefix output/0+_, $(addsuffix .root, $(shell awk -vORS=" " 'NR>1 {print}' sim/temp/converted_widths.txt)))
	@ echo "I'm just a dummy dependency to make Makefile work properly." > $@

# cleans all files created by this section
clean/sim_widths: 
	for file in $(shell awk -vORS=" " 'NR>1 {print}' sim/temp/converted_widths.txt); do $(MAKE) clean/0+_$$file; done
	rm -f sim/temp/*
