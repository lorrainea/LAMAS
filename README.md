# LAMAS: Location-Aware MAtching Statistics

This directory contains the source code, binaries, datasets, and benchmarking scripts for evaluating various implementations of the **LAMAS** algorithm alongside baseline methods (`q-gram`, `ED`, `MS`).

---

## Install LAMAS

cd <algorithm_directory>
./pre-install.sh
make

## Run LAMAS

./lamas -x ../x -y ../y -d 10 -o out

-x <seq_x_filename> input sequence file.
-y <seq_y_filename> input sequence file.
-d <d_range> range of d.
-o <output_filename> output file

