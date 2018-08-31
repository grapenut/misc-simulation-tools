#!/usr/bin/gnuplot

set terminal png

set size square
unset key

set logs xy
set log cb

set pm3d map
set palette model RGB

set xlabel 'Temperature'
set ylabel 'Cooling'
set cblabel 'Metallicity'

set output 'cooling.png'
splot 'dens1.txt' u 3:(($5 - $6)):2 w lines palette

set cblabel 'Density'

set output 'cooling0.png'
splot 'metal0.txt' u 3:(($5 - $6)):1 w lines palette

set output 'cooling1.png'
splot 'metal1.txt' u 3:(($5 - $6)):1 w lines palette

#set ylabel 'Temperature'
#set output 'photoeq.png'
#splot 'output.photoeq.txt' u 1:3:2 w lines palette

#set output 'eden.png'
#plot 'output.photoeq.txt' u 1:4 w lines

#set output 'neutral.png'
#plot 'output.photoeq.txt' u 1:($5/($5+$6)) w lines

