# Variables and paths are supplied by plot_peak.py.
set encoding utf8
set terminal pdfcairo enhanced color font "Helvetica,10" size 3.4in,4.4337in
set output OUTPUT
array bottom[5] = [0.79332709, 0.61962044, 0.44591379, 0.27220714, 0.09850050]
array top[5] = [0.94141184, 0.76770519, 0.59399854, 0.42029189, 0.24658525]
set border lw 1
set tics in scale 1.2,0.6
set xtics mirror font "Helvetica,9"
set ytics mirror font "Helvetica,9"
set mxtics 5
unset mytics
set multiplot
do for [row=1:5] {
    set lmargin at screen 0.16
    set rmargin at screen 0.97
    set bmargin at screen bottom[row]
    set tmargin at screen top[row]
    set xrange [XMIN:XMAX]
    set yrange [-YMAX:YMAX]
    set xtics 20
    set ytics YTICK
    set format y '%g'
    unset xlabel
    set format x ""
    if (row == 5) {
        set format x "%g"
        set xlabel 'r [h^{-1}Mpc]' offset 0,-0.2
    }
    set label 1 sprintf('z = %g', REDSHIFT[row]) at graph 0.05,0.84 left front
    if (row == 1) {
        set title 'ξ(r,z)/D^{2}(z)'
        set key at graph 0.98,0.14 right bottom opaque nobox font "Helvetica,8" \
            maxcols 3 maxrows 1 samplen 1 spacing 0.8
    } else {
        unset title
        unset key
    }
    set arrow 1 from graph 0, first 0 to graph 1, first 0 nohead dt 2 lc rgb "#999999" lw 1 back
    file = sprintf('%s/correlation_%d.dat', TABLES, row)
    plot file using 1:3 with lines lc rgb 'blue' lw 1.0 title 'model', \
         file using 1:4 with lines dt 2 lc rgb 'red' lw 1.0 title 'linear', \
         file using 1:2 with points pt 6 ps 0.35 lw 1 lc rgb 'black' title 'simulation'
    unset arrow 1
    unset label 1
}
unset multiplot
set output
