set terminal png
set output "infected.png"

set xlabel "Time (d)"
set ylabel "Infected ratio (%)"
set grid x y
set ytics nomirror
set tics out

set key left top

# Ajustar rangos si quieres forzar porcentajes
# set y2range [0:100]

set format x "%.0f"
set xtics 10
set logscale y
# set xrange [0: 10]
day = 86400.0

# Si el fichero tiene una cabecera de 1 línea, usamos every ::1 para saltarla.
plot "output/stats-I=   0.50000-R=   0.50000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=0.5", \
"output/stats-I=   0.80000-R=   0.80000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=0.8", \
"output/stats-I=   0.90000-R=   0.90000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=0.9", \
"output/stats-I=   1.00000-R=   1.00000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=1.0", \
"output/stats-I=   1.10000-R=   1.10000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=1.1", \
"output/stats-I=   1.20000-R=   1.20000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=1.2", \
"output/stats-I=   1.50000-R=   1.50000.dat" using ($1/day):2 axes x1y1 with lines lw 2 title "I=R=1.5"
