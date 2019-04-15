#!/bin/bash

./compile.sh

# Handle arguments
inp="HRCornell"
sigma="0.1"
thresh="0.0000001"

if [ $# -ge 1 ] ; then inp=$1; fi
if [ $# -ge 2 ] ; then sigma=$2; fi
if [ $# -ge 3 ] ; then thresh=$3; fi

out=$inp
if [ $# -ge 4 ] ; then out=$4; fi


echo "Input image: $inp"
echo "Sigma: $sigma, Threshold: $thresh"

# Use 'low', 'med' or 'high'
type="med"

# Run 
time ./GSID $inp/$type.ppm "$out/$type-GSID.ppm" $sigma $thresh $inp/normalMap.ppm $inp/depthMap.ppm $inp/shadowMap.ppm $inp/matMap.ppm $inp/colMap.ppm 
