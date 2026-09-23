set encoding utf8
set terminal pdfcairo enhanced color font "Helvetica,8" size 3.4in,1.84in
set output OUTPUT
set logscale xy
set xrange [AMIN:AMAX]
set yrange [1e-7:60]
set border lw 1
set tics in scale 1.2,0.6
set xtics mirror
set ytics mirror
set mxtics 10
set mytics 10
set format x "%g"
set ytics ('10^{-7}' 1e-7, '10^{-5}' 1e-5, '10^{-3}' 1e-3, '10^{-1}' 1e-1, '10^{1}' 10)
set xlabel '{/:Italic a}({/:Italic t})' font 'Helvetica,10' offset 0,0.1
set label 1 '{/:Italic κ}_{/:Italic n}' at screen 0.020,0.54 center rotate by 90 font 'Helvetica,10'
unset key
omit_extra(a) = (abs(a-0.1)<1e-10 || abs(a-0.3)<1e-10 || abs(a-0.6)<1e-10 || abs(a-1)<1e-10) ? 1/0 : a
file = TABLES.'/cumulants.dat'
set multiplot layout 1,2 margins 0.14,0.98,0.22,0.93 spacing 0.05,0
set label 2 sprintf('{/:Italic q} = %g h^{-1}Mpc', Q[1]) at graph 0.10,0.90 left front
plot file using (omit_extra($1)):(abs($2)) with points pt 6 ps 0.42 lw 1 lc rgb 'green', \
     file using (omit_extra($1)):(abs($3)) with points pt 1 ps 0.48 lw 1 lc rgb 'blue', \
     file using (omit_extra($1)):(abs($4)) with points pt 8 ps 0.48 lw 1 lc rgb 'red'
set label 2 sprintf('{/:Italic q} = %g h^{-1}Mpc', Q[4]) at graph 0.10,0.90 left front
set label 11 '{/:Italic κ}_1' at graph 0.80,0.30 right front
set label 12 '{/:Italic κ}_2' at graph 0.80,0.20 right front
set label 13 '{/:Italic κ}_3' at graph 0.80,0.10 right front
set label 21 '' at graph 0.89,0.30 point pt 6 ps 0.42 lw 1 lc rgb 'green' front
set label 22 '' at graph 0.89,0.20 point pt 1 ps 0.48 lw 1 lc rgb 'blue' front
set label 23 '' at graph 0.89,0.10 point pt 8 ps 0.48 lw 1 lc rgb 'red' front
set ytics ('' 1e-7, '' 1e-5, '' 1e-3, '' 1e-1, '' 10)
plot file using (omit_extra($1)):(abs($5)) with points pt 6 ps 0.42 lw 1 lc rgb 'green', \
     file using (omit_extra($1)):(abs($6)) with points pt 1 ps 0.48 lw 1 lc rgb 'blue', \
     file using (omit_extra($1)):(abs($7)) with points pt 8 ps 0.48 lw 1 lc rgb 'red'
unset multiplot
set output
