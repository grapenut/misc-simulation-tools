#!/usr/bin/gnuplot

set terminal png

unset key
set output 't.png'
set yrange [0:6]
set xlabel 'log10 Ionization parameter'
set ylabel 'log10 Temperature [K]'
plot 'data.txt' u 1:(log10($2)) w l 

set output 'h.png'
set key
set yrange [-10:1]
set xlabel 'log10 Ionization parameter'
set ylabel 'log10 n(X)/n(H)'
plot 'data.txt' u 1:3 w l t 'H0', 'data.txt' u 1:4 w l t 'H+', 'data.txt' u 1:5 w l t 'H2'

set output 'he.png'
set ylabel 'log10 n(X)/n(He)'
plot 'data.txt' u 1:6 w l t 'He0', 'data.txt' u 1:7 w l t 'He+', 'data.txt' u 1:8 w l t 'He++'


