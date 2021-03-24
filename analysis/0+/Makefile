# AUSAlib executables
sorter := ~/tools/sorter/build/Sorter
calibrator := ~/tools/calibrator/build/Calibrator
identifier := ~/tools/identifier/build/Identifier

.DEFAULT_GOAL := all
.SECONDARY:

.PHONY: clean cleandat

# removes everything, leaving only the raw data files in data/
clean:
	rm -f sim/raw/*
	rm -f sim/data_i/*
	rm -f sim/data/*
	rm -f match/*
	rm -f id/*
	rm -f output/*
	rm -f -r figures/*

# removes all files associated with %, except for the raw data in data/
# note: also removes any simulated data named %.root
clean/%: 
	@rm -f sim/raw/$(*F).dat
	@rm -f sim/data_i/$(*F).root
	@rm -f sim/data/$(*F).root
	@rm -f match/$(*F)m.root
	@rm -f id/$(*F)mi.root
	@rm -f output/$(*F).root
	@rm -f -r figures/$(*F)
#	@rm -f data/$(*F).root

N = 10

#######################
### tdc_calibration ###
#######################

### temp targets for debugging ###
debug_match: calibrate
	calibrate/mul_compare $(addprefix match/,$(shell ls match))

test: calibrate
	analysis/calibrate_reconstructed output/2455_000.root output/2455_001.root output/2455_002.root

test_all: calibrate $(addprefix output/,$(shell ls data))
	analysis/calibrate_reconstructed $(addprefix output/,$(shell ls data))
### --- ###
		
tdc_calibration: calibrate $(addprefix output/,$(shell ls data))
	analysis/calibrate_true $(addprefix output/,$(shell ls data))	
	analysis/calibrate_reconstructed $(addprefix output/,$(shell ls data))
	hadd -f output/all_events.root output/true_events.root output/reconstructed_events.root

calibrate: analysis/tdc_calibration_reconstructed analysis/tdc_calibration_true ../calibrate/*
	cd analysis/; make calibrate_true; make calibrate_reconstructed
	
#######################
### real root files ###
#######################
### calibration ###
setup/%.cal :
	$(calibrator) -s setup/setup.json -k setup/source.json -d $* caldata/$*/*.root -o setup/ --eloss
calibration : setup/SU.cal setup/SD.cal setup/Det1.cal setup/Det2.cal

### matching ###
match/%m.root : data/%.root matcher.json setup/*
	$(sorter) $< -m matcher.json

# this target allows us to use simulated data from this point on
match/%m.root : sim/data/%.root matcher.json setup/*
	$(sorter) $< -m matcher.json

### identification ###
collectors = --collectors SU:EA:TB,SD:EA:TB,Det1:EA:TF,Det2:EA:TF
target = -t setup/target.json
ions = -i "a"
margin = -m 500keV
id/%i.root : match/%.root
	$(identifier) $(target) $(ions) $(margin) $(collectors) $^ -o id/$*i.root

### analysis ###
analysis/EventAnalysis: analysis/StandardAnalyzer.cpp analysis/main.cpp
	make -C analysis/ EventAnalysis

output/%.root : id/%mi.root analysis/EventAnalysis
	./analysis/EventAnalysis output/$*.root id/$*mi.root match/$*m.root

# Analysis with interference
output/%.root : id/%mi.root analysis/EventAnalysis sim/data_i/%.root
	./analysis/EventAnalysis output/$*.root id/$*mi.root match/$*m.root
	echo sim/data_i/$*.root output/$*.root |root -l ../merger.C

analysis: $(addprefix output/,$(shell ls data))

###############
### figures ###
###############
# the following line allows us to use automatic variables for prerequisites from this point on
.SECONDEXPANSION:

# from here on I'll make heavy use of a few somewhat exotic automatic variables.
# example: pattern: figures/%.pdf, matched by figures/abc/xyz.pdf
# @D --> figures/abc
# *D --> abc
# @F --> xyz.pdf
# *F --> xyz

# compiles a given figure program
../figure_scripts/%: ../figure_scripts/%.cpp
	make -C ../figure_scripts/ $*

# generates a single pdf file for the data source %
figures/%.pdf: output/$$(*D).root ../figure_scripts/$$(*F)
	@ echo "Building $(*F).pdf..."
	@ mkdir -p figures/$(*D)
	@ ../figure_scripts/$(*F) $< $(shell pwd)/figures/$(*D)/

# contains a list of all possible pdf files
figure_pdf = $(addsuffix .pdf,$(basename $(shell find ../figure_scripts/ -maxdepth 1 -name "*.cpp" -printf "%P ")))

temp: 
	@ echo $(figure_pdf)

# concatenates all raw data in data/ and makes all combined figures for them
figs: $(addprefix figures/data/,$(figure_pdf))

# fig/X generates all figures for X. Either data/X.root, ../sim/X.simX, ../sim/X.bahamut must exist
# this command will simulate and/or perform all necessary analysis of the data, and is thus a nice
# point of entry to this whole project 
fig/% : $(addprefix figures/%/,$(figure_pdf)) output/%.root
	@ echo "Built figures/"$*

###################
### simulations ###
###################
# generate the simX configuration file required for simX
sim/%.simX: sim/raw/$$(*F).dat
	@echo "beam: p" > $@
	@echo "target: B11" >> $@
	@echo "-> {" >> $@
	@echo "    g" >> $@
	@echo "    C12 Ex: 17.76MeV" >> $@
	@echo "    -> {" >> $@
	@echo "        generator: PREGEN(file=\""$(shell pwd)/$<"\" unit=\"keV\")" >> $@
	@echo "        a" >> $@
	@echo "        a" >> $@
	@echo "        a" >> $@
	@echo "    }" >> $@
	@echo "}" >> $@

# if the simX file already exists in ../sim/, copy it to sim/
#sim/%.simX: ../sim/%.simX
#	cp $< $@

# generate the raw event data from sim3a
sim/raw/%.dat: ../sim/%.balamuth
	../../sim3a/build/balamuth $(shell cat $<) -N 1000000 -o $(shell pwd)/sim/raw/$(*F).dat	

# simulate the detector setup with simX, using the generated events from sim3a
sim/data/%.root: sim/$$*.simX sim/raw/%.dat
	../sim/simX_fork -c sim/config.json -r sim/$*.simX -o sim/data/$*.root -N 1000000
	@rm $<

# simulate everything solely with simX
sim/data/%.root : ../sim/%.simX
	../sim/simX_fork -c sim/config.json -r $< -o sim/data/$(*F).root -N 1000000
#	@rm $<

###########
# AUTOGEN #
###########
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

###########
# scripts #
###########

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
