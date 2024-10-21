set terminal png size 800,600
set output 'tempi_grafico.png'
set title 'Confronto Tempi Seriali vs Paralleli'
set xlabel 'Dimensione Array'
set ylabel 'Tempo (s)'
set grid
plot 'tempi.csv' using 1:2 with lines title 'Seriale', 'tempi.csv' using 1:3 with lines title 'Parallelo'
