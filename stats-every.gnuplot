set terminal pngcairo size 1200, 1200
infection_rate = 0.08
recovery_rate = 0.01
seed = 42070
set output "stats-every.png"

set grid x y
set ytics
set tics out

unset key

# Ajustar rangos si quieres forzar porcentajes
# set y2range [0:100]

set format x "%.1f"
set xtics 0.1
set ytics 20
set xrange [0:1]
day = 86400.0
file = sprintf('output/stats-I=%10.5f-R=%10.5f-S=%5d.dat', infection_rate, recovery_rate, seed)
set ylabel "Infected ratio (%)"

set multiplot layout 5, 1 \
     title sprintf("I=%f, R=%f", infection_rate, recovery_rate)

set title "Every 1"
# Si el fichero tiene una cabecera de 1 línea, usamos every ::1 para saltarla.
plot file using ($1/day):($4*100) with lines lw 2
set title "Every 10"
plot file every 10 using ($1/day):($4*100) with lines lw 2
set title "Every 100"
plot file every 100 using ($1/day):($4*100) with lines lw 2
set title "Every 1000"
plot file every 1000 using ($1/day):($4*100) with lines lw 2
set title "Every 10000"
set xlabel "Time (d)"
plot file every 10000 using ($1/day):($4*100) with lines lw 2
unset multiplot