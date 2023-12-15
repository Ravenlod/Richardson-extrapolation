#!/bin/bash
echo "RE-exp solver"
SOURCE_DIR="src"
INCLUDE_DIR="include"
OUTPUT_DIR="opt"
OUTPUT_FILE="$OUTPUT_DIR/main_NN"

mkdir $OUTPUT_DIR
if gcc -o "$OUTPUT_FILE" -I "$INCLUDE_DIR" "$SOURCE_DIR"/main_NN.c "$SOURCE_DIR"/odu_NN.c "$SOURCE_DIR"/func_NN.c -lm; then
    "$OUTPUT_FILE"
else
    echo "Compilation failed"
fi
