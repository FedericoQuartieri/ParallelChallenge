set terminal png size 800,600
set output 'tempi_cutoff_max_size.png'
set title 'Confronto Tempi Seriali vs Paralleli'
set xlabel 'Cutoff'
set ylabel 'Tempo (s)'
set grid
set xrange [0:*]
set yrange [0:*]
set datafile separator ','
plot 'tempo_cutoff.csv' using 1:2 with lines title 'Seriale', 'tempo_cutoff.csv' using 1:3 with lines title 'Parallelo'
