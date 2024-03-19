#!/bin/bash

# Define parameters
Gtau_file=./Gtau.dat
Gtau_psd_file=./Gtau_psd.dat
max_it=10 # each iteration performs four projections: density, Hankel, PSD, PSD (top submatrix)
flavors=2 
ntau=2001 # must be odd
tol=1e-12 
verbose=0 # 0 or 1
app=../build/dyproject

# Run the program with the specified parameters
$app --Gtau "$Gtau_file" --Gtau_psd "$Gtau_psd_file" --max_it "$max_it" --flavors "$flavors" --ntau "$ntau" --tol "$tol" --verbose "$verbose"

