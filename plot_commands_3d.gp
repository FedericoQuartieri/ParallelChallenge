set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'plot_3d_cutoff_size.png'
set title '3D Comparison of Parallel Times (Cutoff and Array Size)'
set xlabel 'Cutoff'
set ylabel 'Array Size'
set zlabel 'Parallel Time (s)'
set grid
set datafile separator ','
set xrange [0:7]
set yrange [100000:1e+06]
set zrange [0:*]
set pm3d at s
set dgrid3d 100,100
splot 'tempi_3d_interpolated.csv' using 1:2:3 with pm3d title 'Parallel Time'
