#!/usr/bin/gnuplot

set terminal png

set size square
unset key

set logs xy
set log cb

set pm3d map
set palette model RGB

set xlabel 'Ionization parameter'
set ylabel 'Thermal Time (years)'
set cblabel 'Metallicity'

set output 'thermal.png'
splot 'output.photoeq.txt' u 1:($13/3e7):2 w lines palette

set ylabel 'Temperature'
set output 'photoeq.png'
splot 'output.photoeq.txt' u 1:3:2 w lines palette

#set output 'eden.png'
#plot 'output.photoeq.txt' u 1:4 w lines

#set output 'neutral.png'
#plot 'output.photoeq.txt' u 1:($5/($5+$6)) w lines

