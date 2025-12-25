set terminal png
set output "infected-processed.png"
set datafile separator ','
set key autotitle columnhead
set xlabel "Time (d)"
set ylabel "Infected ratio (%)"
set grid x y
set ytics nomirror
set tics out

set key right top

# Ajustar rangos si quieres forzar porcentajes
# set y2range [0:100]

set format x "%.0f"
set xtics 10
# set logscale y
# set xrange [0: 10]
day = 86400.0

# Si el fichero tiene una cabecera de 1 línea, usamos every ::1 para saltarla.
plot "processed/I=R=0.7.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=0.7", \
"processed/I=R=0.8.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=0.8", \
"processed/I=R=0.9.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=0.9", \
"processed/I=R=1.0.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=1.0", \
"processed/I=R=1.1.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=1.1", \
"processed/I=R=1.2.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=1.2", \
"processed/I=R=1.3.csv" using ($1/day):2:3 axes x1y1 with lines lw 2 title "I=R=1.3"
