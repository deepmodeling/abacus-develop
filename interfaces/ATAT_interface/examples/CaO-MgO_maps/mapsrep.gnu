set terminal pngcairo enhanced size 1200,900 font "Arial,12"
set output 'mapsrep.png'

set xlabel "concentration"
set ylabel "energy"
set title "Fitted Energies"
plot [0:1] \
'predstr.out' u 1:3 t "predicted" w p pt 2, \
'fit.out' u 1:3 t  "known str" w p pt 1, \
'gs.out' u 1:3 t "known gs" w linesp lt 3
pause -1
set title "Calculated Energies"
plot [0:1] \
'fit.out' u 1:2 t  "known str" w p pt 5, \
'gs.out' u 1:2 t "known gs" w linesp lt 6
pause -1
set title "Calculated and Fitted Energies"
plot [0:1] \
'fit.out' u 1:2 t  "calculated" w p pt 5, \
'fit.out' u 1:3 t  "fitted" w p pt 1 
pause -1
set xlabel "diameter"
set ylabel "energy"
set title "ECI vs cluster diameter"
set nokey
set xzeroaxis
set xtics ("pair" 0,"5" 5,"10" 10,"15" 15, "trip" 20,"5" 25,"10" 30,"15" 35, "quad" 40,"5" 45,"10" 50,"15" 55)
plot [-5:60] 'clusinfo.out' u (($1-2)*20.+$2):($4)
pause -1
set title "Residuals of the fit (same order as in fit.out)"
set xtics autofreq
set ylabel "energy"
set xlabel "line number in fit.out"
plot 'fit.out' u 4
pause -1
set title "Fitted vs Calculated Energies"
set ylabel "predicted energy"
set xlabel "actual energy"
set nokey
plot \
'fit.out' u 2:3 w p pt 1,x
pause -1
