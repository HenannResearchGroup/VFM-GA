#!/bin/bash
module load intel
cd \GA_fort
ifort constants.f90 linalg.f90 fem_basics.f90 fem_const.f90 fem_compute_el.f90 stability.f90 fem_compute_global.f90 fem_preprocess.f90 print_write.f90 genetic_algo.f90 main.f90 -o run_ga.exe