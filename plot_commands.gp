set terminal png size 800,600
set output 'tempi_size_11.png'
set title 'Comparison of Serial vs Parallel Times with cutoff 11
set xlabel 'Array Size'
set ylabel 'Time (s)'
set grid
set xrange [100000:*]
set yrange [0:*]
set datafile separator ','
plot 'tempi_size_11.csv' using 1:2 with lines title 'Serial', 'tempi_size_11.csv' using 1:3 with lines title 'Parallel'
