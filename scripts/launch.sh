#!bin/bash

echo "RE-exp solver"

out=$(gcc -o var/main_NN -I include/ src/main_NN.c src/odu_NN.c src/func_NN.c -lm)

echo $(out)