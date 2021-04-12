# AUSAlib executables. If globally installed, replace with "./"
sorter := ~/tools/sorter/build/Sorter
calibrator := ~/tools/calibrator/build/Calibrator
identifier := ~/tools/identifier/build/Identifier

.DEFAULT_GOAL := help
.SECONDARY:
.PHONY: clean clean_build

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
clean: clean_build
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
balamuth_i := ../../sim3a_i/build/balamuth

# if any of the dependent executables doesn't exist, we compile them
../../simX_fork/build/simX_fork/build/simX:
	make -C ../../simX_fork/build/
	
../../sim3a/build/balamuth: ../../sim3a/balamuth/*
	make -C ../../sim3a/build/

../../sim3a_i/build/balamuth: ../../sim3a_i/balamuth/*
	make -C ../../sim3a_i/build/

# removes all executables
clean_build: 
	@ rm -f -r ../figure_scripts/build/*
	@ rm -f -r ../special_figure_scripts/build/*
	@ rm -f -r ../calibrate/build/*
	@ rm -f -r analysis/build/*

#######################
### tdc_calibration ###
#######################
calibrate_rec := analysis/build/tdc_calibration_reconstructed
calibrate_true := analysis/build/tdc_calibration_true

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

output/%.root : match/%m.root id/%mi.root $(event_analysis)
	$(event_analysis) output/$*.root id/$*mi.root match/$*m.root

# Analysis with interference
output/%.root : id/%mi.root $(event_analysis) sim/data_i/%.root
	$(event_analysis) output/$*.root id/$*mi.root match/$*m.root
	echo sim/data_i/$*.root output/$*.root |root -l ../merger.C

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
../scripts/special/build/%: ../scripts/special/%.cpp
	make -C ../scripts/special/build/ $*

# creates an interactive 3D Dalitz plot of %. mainly used for debug
Dalitz_3D/%: output/%.root ../scripts/special/build/Dalitz_3D
	@ echo "Creating interactive 3D Dalitz plot..."
	../scripts/special/build/Dalitz_3D $<

# creates two histograms, one for % and another for output/all_events.
# it then binwise divides each histogram to give the output histogram, showing the deviations
figures/%/sim_compare_raw.pdf: output/%.root output/all_events.root ../scripts/special/build/sim_compare
	@ echo "Building "$@"..."
	../scripts/special/build/sim_compare output/all_events.root $< $(shell pwd)/figures/$*/

# same as above, but now with a cut on bin values to enhance smaller details and fix the z axis
figures/%/sim_compare_cut.pdf: output/%.root output/all_events.root ../scripts/special/build/sim_compare
	@ echo "Building "$@"..."
	../scripts/special/build/sim_compare output/all_events.root $< $(shell pwd)/figures/$*/ 3

# this is just for easy access to the two commands above. simply use "make sim_compare/X" to create both sim_compare plots for output/X.root
sim_compare/%: figures/%/sim_compare_raw.pdf figures/%/sim_compare_cut.pdf;

### thesis specific figures ###
# plots before / after of the TDC values
tdc_plots: ../scripts/special/build/tdc_plot $(addprefix output/,$(shell ls data))
	@ echo "Creating figures/article/tdc_plot_before.pdf and figures/article/tdc_plot_after.pdf"
	$< $(addprefix output/,$(shell ls data)) $(shell pwd)/figures/article/

# illustrates how the peaks depend on the multiplicity
mul_compare: ../scripts/special/build/mul_compare $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data))))
	@ echo "Creating figures/article/mul_compare.pdf"
	$< $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data)))) $(shell pwd)/figures/article/

mul_compare_debug: ../scripts/special/build/mul_compare $(addsuffix m.root, $(addprefix match/, $(basename $(shell ls data))))
	@ echo "Creating figures/article/mul_compare.pdf"
	$< match/2455_000m.root $(shell pwd)/figures/article/

###################
### SIMULATIONS ###
###################
# number of events generated by all simulations
N := 1000000

# generate the simX configuration file required for simX
# this is generated on the fly whenever it is needed, and is deleted afterwards
sim/%.simX: sim/raw/$$(*F).dat
	@echo "beam: p" > $@
	@echo "target: B11" >> $@
	@echo "-> {" >> $@
	@echo "    generator: PREGEN(file=\""$(shell pwd)/$<"\" unit=\"keV\")" >> $@
	@echo "    a" >> $@
	@echo "    a" >> $@
	@echo "    a" >> $@
	@echo "}" >> $@

# generate the raw event data from sim3a
sim/raw/%.dat: sim/configs/%.balamuth $(balamuth)
	$(balamuth) $(shell cat $<) -N $(N) -o $(shell pwd)/sim/raw/$(*F).dat	

# simulate the detector setup with simX, using the generated events from sim3a
sim/data/%.root: sim/$$*.simX sim/raw/%.dat $(simX)
	$(simX) -c sim/config.json -r sim/$*.simX -o sim/data/$*.root -N $(N)
	@ rm $<

# simulate everything solely with simX
sim/data/%.root : sim/configs/%.simX $(simX)
	$(simX) -c sim/config.json -r $< -o sim/data/$(*F).root -N $(N)

###############
### AUTOGEN ###
###############
directory = ~/projects/B11/autogen/$(notdir $(shell pwd))
autogen :
	+make autofig/auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC)
	rm $(directory)/sim/data_i/auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC).*.root
	rm $(directory)/sim/data/auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC).*.root
	rm $(directory)/match/auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC).*.root
	rm $(directory)/id/auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC).*.root
	rm $(directory)/output/auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC).*.root

$(directory)/sim/data_i/auto_%.root :
	cd ../sim; ./balamuth_i --JC $(JC) --JBe 2+ -L 3 --ExC 16.6MeV --L2 1 --phase 6 --ratio 0.7 -N 1000000 -o $(directory)/sim/data_i/auto_'$*'.root --seed 50$(subst .,,$(suffix $*)) --coulomb-radius $(coulomb-radius)fm --r0Be $(r0Be)fm --gBe $(gBe)'keV^(1/2)' --ExBe $(ExBe)keV
	#cd ../sim; ./balamuth_i --JC $(JC) --JBe 2+ -L 3 --ExC 16.6MeV --L2 1 --phase 6 --ratio 0.7 -N 1000000 -o $(directory)/sim/data_i/auto_'$*'.root --seed 50$(subst .,,$(suffix $*)) --r0Be $(r0Be)fm --gBe $(gBe)'keV^(1/2)' --ExBe $(ExBe)keV

$(directory)/sim/data/%.root : $(directory)/sim/data_i/%.root
	cd ../sim; python3 simX_factory.py $(directory)/sim/data_i/$*.root > $(directory)/sim/data/$*.simX
	cd ../sim; ./simX_fork -c $(shell pwd)/sim/config.json -r $(directory)/sim/data/$*.simX -o $(directory)/sim/data/$*.root -N 1000000 -s $*
	rm $(directory)/sim/data/$*.simX

$(directory)/match/%m.root : $(directory)/sim/data/%.root matcher.json setup/*
	$(sorter) $< -m matcher.json -o $(directory)/match/$*m.root

$(directory)/id/%i.root : $(directory)/match/%.root
	$(identifier) $(target) $(ions) $(margin) $(collectors) $^ -o $(directory)/id/$*i.root

$(directory)/output/%.root : $(directory)/id/%mi.root analysis/EventAnalysis $(directory)/sim/data_i/%.root
	./analysis/EventAnalysis $(directory)/output/$*.root $(directory)/id/$*mi.root $(directory)/match/$*m.root
	echo $(directory)/sim/data_i/$*.root $(directory)/output/$*.root |root -l ../merger.C

$(directory)/figures/%.pdf : $(addprefix $(directory)/,$$(call fun1,%)) ../figure_scripts/$$(call fun,%).py $(addprefix $(directory)/,$(call simfiles,auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC),$(N)))
	mkdir -p $(directory)/figures/$(subst /,,$(dir $*))
	echo $(addprefix $(directory)/,$(call fun1,$*)) | python3 ../figure_scripts/$(notdir $(basename $@)).py $(directory)/figures/$(subst /,,$(dir $*))

autofig/% : $(addprefix $(directory)/,$(addprefix figures/%/,$(figure_pdf)))
	echo $(addprefix $(directory)/,$(addprefix figures/$*/,$(figure_pdf)))

automaker :
	cat ../list2 |python3 ../automaker.py $(shell pwd)

###############
### SCRIPTS ###
###############

test.pdf : ../scripts/test.py figures/data_tdc_cut/dalitz_forfit.pdf
	python3 ../scripts/test.py $(shell pwd) > test.log
test_automake.pdf : ../scripts/test_automake.py figures/data_tdc_cut/dalitz_forfit.pdf
	python3 ../scripts/test_automake.py $(shell pwd) > test_automake.log

surf.pdf : ../scripts/surf.py
	python3 ../scripts/surf.py $(shell pwd)

test_noecut.pdf : ../scripts/test_noecut.py
	python3 ../scripts/test_noecut.py $(shell pwd) > test_noecut.log
	
all.pdf : ../scripts/all.py
	python3 ../scripts/all.py $(shell pwd) > all.log


%.pdf : ../scripts/%.py
	python3 ../scripts/$*.py $(shell pwd)
