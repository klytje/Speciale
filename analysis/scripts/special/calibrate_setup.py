import sys

n_in = sys.argv[1]
n_out = sys.argv[2]

f_in = open(n_in, "r")
f_out = open(n_out, "w")

lines = f_in.readlines()

for line in lines:
    if "calibration" in line: 
        continue
    f_out.write(line)
