# AUSAlib executables
sorter := ~/tools/sorter/build/Sorter
calibrator := ~/tools/calibrator/build/Calibrator
identifier := ~/tools/identifier/build/Identifier

.DEFAULT_GOAL := all
.SECONDARY:

.PHONY: clean cleandat tdc_calibration calibrate mul_compare
clean :
	rm -f sim/data_i/*
	rm -f sim/data/*
	rm -f match/*
	rm -f id/*

cleandat :
	rm -f sim/data_i/*.dat

N = 10

#######################
### tdc_calibration ###
#######################
debug_match: ../Calibrate/mul_compare
	cd ../Calibrate; make mul_compare
	../Calibrate/mul_compare $(addprefix match/,$(shell ls match))
		
tdc_calibration: calibrate $(addprefix output/,$(shell ls data))
	analysis/calibrate_true $(addprefix output/,$(shell ls data))	
	analysis/calibrate_reconstructed $(addprefix output/,$(shell ls data))
	hadd -f output/all_events.root $^

calibrate:
	cd analysis/; make calibrate_true; make calibrate_reconstructed
	
test: calibrate
	analysis/calibrate_reconstructed output/2414_000.root output/2414_001.root output/2414_002.root

test_all: calibrate $(addprefix output/,$(shell ls data))
	analysis/calibrate_reconstructed $(addprefix output/,$(shell ls data))

merged/%.root : ../tdcmerger.C match/%m.root output/%.root
	echo $* |root -b ../tdcmerger.C

all_merged : $(addprefix merged/,$(shell ls data))

# Calibrator
setup/%.cal :
	$(calibrator) -s setup/setup.json -k setup/source.json -d $* caldata/$*/*.root -o setup/ --eloss
calibration : setup/SU.cal setup/SD.cal setup/Det1.cal setup/Det2.cal

# Analysis code
analysis/EventAnalysis : analysis/main.cpp analysis/StandardAnalyzer.cpp
	cd analysis; make EventAnalysis

# Sorting / Matching
match/%m.root : sim/data/%.root matcher.json setup/*
	$(sorter) $< -m matcher.json
match/%m.root : data/%.root matcher.json setup/*
	$(sorter) $< -m matcher.json

# Identification
collectors = --collectors SU:EA:TB,SD:EA:TB,Det1:EA:TF,Det2:Ea:TF
target = -t setup/target.json
ions = -i "a"
margin = -m 500keV
id/%i.root : match/%.root
	$(identifier) $(target) $(ions) $(margin) $(collectors) $^ -o id/$*i.root

# Analysis with interference
output/%.root : id/%mi.root analysis/EventAnalysis sim/data_i/%.root
	./analysis/EventAnalysis output/$*.root id/$*mi.root match/$*m.root
	echo sim/data_i/$*.root output/$*.root |root -l ../merger.C
# Analysis
output/%.root : id/%mi.root analysis/EventAnalysis
	./analysis/EventAnalysis output/$*.root id/$*mi.root match/$*m.root

# Figures
output_data = $(addprefix output/,$(shell ls data/ ))

# example: figures/abc/xyz.pdf
# @D --> figures/abc
# *D --> abc
# @F --> xyz.pdf
# *F --> xyz
figures/%.pdf:
#	ifeq ("$(wildcard output/$(*D))","")
#		$(error Could not find root file: output/$(*D).root)
#	endif
	@ echo "Building $(*F).pdf..."
	@ mkdir -p figures/$(*D)
	@ echo output/$(*D).root | python3 ../figure_scripts/$(*F).py $(shell pwd)/$(@D)

figures/data/%.pdf : $(output_data) ../figure_scripts/%.py
	@ mkdir -p figures/data
	@ echo $(output_data) | python3 ../figure_scripts/$*.py $(shell pwd)/figures/data

figures/data_tdc_cut/%.pdf : output/data_tdc_cut.root ../figure_scripts/%.py
	@ mkdir -p figures/data_tdc_cut
	@ echo output/data_tdc_cut.root | python3 ../figure_scripts/$*.py $(shell pwd)/figures/data_tdc_cut
	
figures/data_tdc_cut_inverse/%.pdf : output/data_tdc_cut.root ../figure_scripts/%.py
	@ mkdir -p figures/data_tdc_cut_inverse
	@ echo output/data_tdc_cut_inverse.root | python3 ../figure_scripts/$*.py $(shell pwd)/figures/data_tdc_cut_inverse

simfiles = $(shell echo output/$(1).{1..$(2)}.root)

fun = $(notdir $(1))
fun1 = $(call simfiles,$(subst /,,$(dir $(1))),$N)
.SECONDEXPANSION:
figures/%.pdf : $$(call fun1,%) ../figure_scripts/$$(call fun,%).py# figures/dataframe.py
	mkdir -p figures/$(subst /,,$(dir $*))
	echo $(call fun1,$*) | python3 ../figure_scripts/$(notdir $(basename $@)).py $(shell pwd)/figures/$(subst /,,$(dir $*))

# All
figure_pdf = $(addsuffix .pdf,$(basename $(shell ls ../figure_scripts/ | grep py | grep -v dataframe | grep -v pycache)))

figs: $(addprefix figures/data/,$(figure_pdf))

fig/% : $(addprefix figures/%/,$(figure_pdf))
	# Built figures/$*

# Simulate
.SECONDEXPANSION:
sim/data/%.root : $$(basename ../sim/%).simX
	cd ../sim; ./simX_fork -c $(shell pwd)/sim/config.json -r $(notdir $<) -o $(shell pwd)/sim/data/$*.root -N 1000000 -s $*

# Simulate scrambled
.SECONDEXPANSION:
sim/data/random.%.root : ../sim/balamuth.template
	cd ../sim; python3 rgen.py > $(shell pwd)/sim/data/random.$*.dat
	cd ../sim; python3 simX_factory.py $(shell pwd)/sim/data/random.$*.dat random.template > $(shell pwd)/sim/data/random.$*.simX
	cd ../sim; ./simX_fork -c $(shell pwd)/sim/config.json -r $(shell pwd)/sim/data/random.$*.simX -o $(shell pwd)/sim/data/random.$*.root -N 1000000 -s $*
	rm $(shell pwd)/sim/data/random.$*.dat $(shell pwd)/sim/data/random.$*.simX

# Simulate balamuth
.SECONDEXPANSION:
sim/data/%.root : $$(basename ../sim/%).balamuth
	cd ../sim; ./balamuth_i $(shell cat $(basename ../sim/$*).balamuth) -N 1000000 -o $(shell pwd)/sim/data/$*.dat --seed $(subst .,,$(suffix $*))
	cd ../sim; python3 simX_factory.py $(shell pwd)/sim/data/$*.dat > $(shell pwd)/sim/data/$*.simX
	cd ../sim; ./simX_fork -c $(shell pwd)/sim/config.json -r $(shell pwd)/sim/data/$*.simX -o $(shell pwd)/sim/data/$*.root -N 1000000 -s $*
	rm $(shell pwd)/sim/data/$*.dat $(shell pwd)/sim/data/$*.simX


.SECONDEXPANSION:
sim/data_i/%.root : $$(basename ../sim/%).balamuth_i
	cd ../sim; ./balamuth_i $(shell cat $(basename ../sim/$*).balamuth_i) -N 1000000 -o $(shell pwd)/sim/data_i/$*.root --seed $(subst .,,$(suffix $*))

.SECONDEXPANSION:
sim/data/%.root : sim/data_i/%.root
	cd ../sim; python3 simX_factory.py $(shell pwd)/sim/data_i/$*.root > $(shell pwd)/sim/data/$*.simX
	cd ../sim; ./simX_fork -c $(shell pwd)/sim/config.json -r $(shell pwd)/sim/data/$*.simX -o $(shell pwd)/sim/data/$*.root -N 1000000 -s $*
	rm $(shell pwd)/sim/data/$*.simX

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

.SECONDEXPANSION:
$(directory)/sim/data_i/auto_%.root :
	cd ../sim; ./balamuth_i --JC $(JC) --JBe 2+ -L 3 --ExC 16.6MeV --L2 1 --phase 6 --ratio 0.7 -N 1000000 -o $(directory)/sim/data_i/auto_'$*'.root --seed 50$(subst .,,$(suffix $*)) --coulomb-radius $(coulomb-radius)fm --r0Be $(r0Be)fm --gBe $(gBe)'keV^(1/2)' --ExBe $(ExBe)keV
	#cd ../sim; ./balamuth_i --JC $(JC) --JBe 2+ -L 3 --ExC 16.6MeV --L2 1 --phase 6 --ratio 0.7 -N 1000000 -o $(directory)/sim/data_i/auto_'$*'.root --seed 50$(subst .,,$(suffix $*)) --r0Be $(r0Be)fm --gBe $(gBe)'keV^(1/2)' --ExBe $(ExBe)keV

.SECONDEXPANSION:
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

.SECONDEXPANSION:
$(directory)/figures/%.pdf : $(addprefix $(directory)/,$$(call fun1,%)) ../figure_scripts/$$(call fun,%).py $(addprefix $(directory)/,$(call simfiles,auto_$(coulomb-radius)_$(r0Be)_$(gBe)_$(ExBe)_$(JC),$(N)))
	mkdir -p $(directory)/figures/$(subst /,,$(dir $*))
	echo $(addprefix $(directory)/,$(call fun1,$*)) | python3 ../figure_scripts/$(notdir $(basename $@)).py $(directory)/figures/$(subst /,,$(dir $*))

.SECONDEXPANSION:
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
