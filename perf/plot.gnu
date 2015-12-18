set terminal png nocrop enhanced font "arial,8" 
set output 'results.png'
set grid nopolar
set grid xtics nomxtics ytics nomytics noztics nomztics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set style data lines
set title "Time to compute solution, function of number of processors and matrix size" 
set xlabel "Nbr of processors"
set ylabel "Matrix size (n)"
set zlabel "Time to compute"
## Last datafile plotted: "$grid"
splot 'res04.dat' using 1:2:3 with points title "Computation time"
