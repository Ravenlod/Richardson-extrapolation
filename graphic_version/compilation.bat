@echo off

set COMPILER_COMMAND=gcc exp.c -o exp

%COMPILER_COMMAND%

if %errorlevel% equ 0 (
    echo Success!
    start "" exp
    call .\gnuplot\bin\gnuplot .\plot_script.gnu
) else (
    echo Error.
)
exit