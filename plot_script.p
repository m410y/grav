set terminal png size 500, 500
set output 'grav.png'
set xrange[0:1]
set yrange[0:1]
set title 'Torus'
plot 'output.txt' using 2:3 with dots
