set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'plot_3d_cutoff_size.png'
set title 'Confronto Tempi Paralleli in 3D (Cutoff e Dimensione Array)'
set xlabel 'Cutoff'
set ylabel 'Dimensione Array'
set zlabel 'Tempo Parallelo (s)'
set grid
set datafile separator ','
set xrange [0:8]
set yrange [0 :1000000]
set zrange [0:*]
set pm3d at s
set dgrid3d 100,100
splot 'tempi_3d.csv' using 1:2:4 with pm3d title 'Parallel Time'
