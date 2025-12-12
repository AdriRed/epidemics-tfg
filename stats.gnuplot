set terminal png
set output "stats.png"

set xlabel "Time (d)"
set ylabel "Rate"
set y2label "Infected ratio (%)"
set grid x y

set ytics nomirror
set y2tics
set tics out

set key center top

# ajustar rangos si quieres forzar porcentajes
#set y2range [0:100]

set format x "%.0f"
set xtics 10

day = 86400.0
file = 'output/stats-I= 1.000000000-R= 0.010000000.dat'
# Si el fichero tiene una cabecera de 1 línea, usamos every ::1 para saltarla.
# Ajusta or quita every ::1 si no hay cabecera.
plot file every ::1 using ($1/day):2 axes x1y1 with lines lw 2 title "Infection Rate", \
     file every ::1 using ($1/day):3 axes x1y1 with lines lw 2 title "Recovery Rate", \
     file every ::1 using ($1/day):($4*100) axes x1y2 with lines lw 2 title "Infected population (%)"