set term wxt 1
set xlabel "x"
set xlabel "y"
plot "result.txt" using 1:2 with linespoints ls 1 title 'y1(x)'
set term wxt 2
set xlabel "x"
set xlabel "y"
plot "result.txt" using 1:3 with linespoints ls 1 title 'y2(x)'
set term wxt 3
set xlabel "x"
set xlabel "y"
plot "result.txt" using 1:4 with linespoints ls 1 title 'y3(x)'
pause -1