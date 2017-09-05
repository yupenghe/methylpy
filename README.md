# methylpy

Welcome to the home page of methylpy a whole genome bisulfite sequencing analysis pipeline.

## Unit test for DMRfind
```
cd methylpy/test/
python2.7 DMRfind_test.py
```

## Compiling rms.cpp
* Most cases
g++ -O3 -l gsl -l gslcblas -o run_rms_tests.out rms.cpp
* Ubuntu 16.04
g++ -o run_rms_tests.out rms.cpp `gsl-config --cflags —libs`
