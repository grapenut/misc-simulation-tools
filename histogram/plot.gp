#!/usr/bin/gnuplot

set terminal postscript eps enhanced solid color "Arial" 20
set output 'histogram.ps'

#set terminal x11 enhanced

set size square

set key top right

set key samplen 2
set key spacing 0.9

set logs

set format y "%T"
set format x "%T"

set border lw 1.0

set tics front scale 1.5
set mxtics 10
set mytics 10

shaded_color = "#DDDDDD"
fluid_color = "black"
sn_colors = "red orange yellow green cyan blue violet"

ymin = 1e-9
set xrange [1e-3:1e5]
set yrange [ymin:2e-3]

set xlabel 'Log Number Density [cm^{-3}]'
set ylabel 'Log Metallicity [M_Z / M_{total}]'

file = 'hist_data_1753.txt'

#plot file u 2:($5-sqrt($6) > 0 ? ($5-sqrt($6)) : ymin):($5+sqrt($6)) w filledcurves t '' lw 3 lc rgb shaded_color, file u 2:5 w line lw 3 lc rgb fluid_color t 'Fluid', file u 2:9 w line lw 3 lc rgb word(sn_colors,1) t 'SN 1', file u 2:11 w line lw 3 lc rgb word(sn_colors,2) t 'SN 2', file u 2:13 w line lw 3 lc rgb word(sn_colors,3) t 'SN 3', file u 2:15 w line lw 3 lc rgb word(sn_colors,4) t 'SN 4', file u 2:17 w line lw 3 lc rgb word(sn_colors,5) t 'SN 5', file u 2:19 w line lw 3 lc rgb word(sn_colors,6) t 'SN 6', file u 2:21 w line lw 3 lc rgb word(sn_colors,7) t 'SN 7'
plot file u 2:5 w line lw 3 lc rgb fluid_color t 'Fluid', file u 2:9 w line lw 3 lc rgb word(sn_colors,1) t 'SN 1', file u 2:11 w line lw 3 lc rgb word(sn_colors,2) t 'SN 2', file u 2:13 w line lw 3 lc rgb word(sn_colors,3) t 'SN 3', file u 2:15 w line lw 3 lc rgb word(sn_colors,4) t 'SN 4', file u 2:17 w line lw 3 lc rgb word(sn_colors,5) t 'SN 5', file u 2:19 w line lw 3 lc rgb word(sn_colors,6) t 'SN 6', file u 2:21 w line lw 3 lc rgb word(sn_colors,7) t 'SN 7'

pause -1

#unset logs
#set log x

set ylabel 'Discrete Error [M_{Z,\rm{fluid}} / M_{Z,{\rm particles})]'
set format y '%.0f%%'

set output 'fluid_particle_error.ps'
set yrange [1:50]

plot file u 2:((($4-$7)/$4)*100.0) w line t ''

pause -1

set output 'fluid_particle_mass.ps'
set format y '%T'
set yrange [0.0001:10]
set ylabel 'Mass [M_\odot]'
plot file u 2:($4/2e33) w line t 'Fluid Metal Mass' lw 3 lc rgb 'black', file u 2:($7/2e33) w line t 'Particle Mass' lw 3 lc rgb 'red'
pause -1





