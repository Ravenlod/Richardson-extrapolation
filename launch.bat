@echo off
echo RE-exp solver

set "SOURCE_DIR=src"
set "INCLUDE_DIR=include"
set "OUTPUT_DIR=opt"
set "OUTPUT_FILE=%OUTPUT_DIR%\main_NN.exe"

mkdir "%OUTPUT_DIR%"
gcc -o "%OUTPUT_FILE%" -I "%INCLUDE_DIR%" "%SOURCE_DIR%\main_NN.c" "%SOURCE_DIR%\odu_NN.c" "%SOURCE_DIR%\func_NN.c" -lm

if %errorlevel% equ 0 (
    echo Compilation successful
    "%OUTPUT_FILE%"
) else (
    echo Compilation failed
)
pause -1