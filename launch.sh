#!/bin/bash

echo "RE-exp solver"

gcc -o opt/main_NN -I include/ src/main_NN.c src/odu_NN.c src/func_NN.c -lm

opt/main_NN