set terminal png
infection_rate = 1.0
recovery_rate = 1.0
set output sprintf("stats-I=%f-R=%f.png", infection_rate, recovery_rate)

set xlabel "Time (d)"
set ylabel "Rate"
set y2label "Infected ratio (%)"
set grid x y
set title sprintf("I=%f, R=%f", infection_rate, recovery_rate)
set ytics nomirror
set y2tics
set tics out

set key center top

# Ajustar rangos si quieres forzar porcentajes
# set y2range [0:100]

set format x "%.0f"
set xtics 1
day = 86400.0
file = sprintf('output/stats-I= %.9f-R= %.9f.dat', infection_rate, recovery_rate)

# Si el fichero tiene una cabecera de 1 línea, usamos every ::1 para saltarla.
plot file every ::1 using ($1/day):2 axes x1y1 with lines lw 2 title "Infection Rate", \
     file every ::1 using ($1/day):3 axes x1y1 with lines lw 2 title "Recovery Rate", \
     file every ::1 using ($1/day):($4*100) axes x1y2 with lines lw 2 title "Infected population (%)"
