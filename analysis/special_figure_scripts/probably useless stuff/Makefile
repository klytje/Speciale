# I used this once, and it may become again relevant later, so here it is. 
# It requires two files to be placed in this very folder: true_aligned_peaks.root and reconstructed_aligned_peaks.root
# they can be created by enabling the lines "save(&data, X)" in analysis/tdc_calibration_reconstructed.cpp and analysis/tdc_calibration_true.cpp
# after they have been enabled, simply run "make tdc_calibration" and they'll appear in the main experiment folder.

# first the FT and BT peaks are repaired and shifted. then they are normalized and placed on top of each other, such that their shapes can be compared. this may help reveal why a given TDC calibration is poor (more specifically, why the sigma values are so large on the gaussian cut figures)

aligned.root: true_aligned_peaks.root reconstructed_aligned_peaks.root
	hadd $@ $^

# the only difference between these two targets is the amount FT and BT are shifted by (different offsets for different experiments)
0+: aligned.root
	root test_0+.C
	
2-: aligned.root
	root test_2-.C
