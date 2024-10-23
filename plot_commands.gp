set terminal png size 800,600
set output 'tempi_cutoff_max_size.png'
set title 'Comparison of Serial vs Parallel Times'
set xlabel 'Cutoff'
set ylabel 'Time (s)'
set grid
set xrange [0:*]
set yrange [0:*]
set datafile separator ','
plot 'tempo_cutoff.csv' using 1:2 with lines title 'Serial', 'tempo_cutoff.csv' using 1:3 with lines title 'Parallel'
