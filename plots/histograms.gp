# Transition histograms and sampled trajectories, rendered with Cairo text.
set encoding utf8
set terminal pdfcairo enhanced color font 'Helvetica,10' size 11.2535in,8.3228in
set output OUTPUT
array colours[4] = ['#ff0000', '#000000', '#0000ff', '#00ff00']
array pale[4] = ['#ffb3b3', '#b3b3b3', '#b3b3ff', '#b3ffb3']
array times[4] = [0.1, 0.3, 0.6, 1.0]
array bottom[4] = [0.72115, 0.49880, 0.27885, 0.05890]
array top[4] = [0.91106, 0.69110, 0.46880, 0.24880]
array left[4] = [0.30422, 0.46322, 0.62222, 0.78122]
array right[4] = [0.45122, 0.61022, 0.76922, 0.92822]
set border lw 1
set tics in scale 1.2,0.6
set xtics mirror
set ytics mirror
unset key
set arrow 100 from screen 0.3961,0.962 to screen 0.8635,0.962 head filled size screen 0.010,15,45 front lw 1
set label 100 'increasing time' at screen 0.6298,0.978 center font 'Helvetica,12'
set arrow 101 from screen 0.955,0.801 to screen 0.955,0.171 head filled size screen 0.010,15,45 front lw 1
set label 101 'decreasing initial separation {/:Italic q}' at screen 0.975,0.486 center rotate by -90 font 'Helvetica,12'
set multiplot
do for [row=1:4] {
    set lmargin at screen 0.063
    set rmargin at screen 0.23022
    set bmargin at screen bottom[row]
    set tmargin at screen top[row]
    unset logscale y
    unset grid
    unset title
    set xrange [0:1]
    set yrange [-40:40]
    set xtics 0.25
    set mxtics 5
    set ytics 20
    set mytics 5
    set format y '%g'
    set format x ''
    unset xlabel
    if (row == 4) {
        set format x '%g'
        set xlabel '{/:Italic a}({/:Italic t})' font 'Helvetica,12' offset 0,-0.2
    }
    set ylabel '{/:Italic r} − {/:Italic q}  [h^{-1}Mpc]' offset -1.0,0
    set label 1 sprintf('{/:Italic q} = %g h^{-1}Mpc', Q[row]) at graph 0.13,0.88 left front tc rgb colours[row] font 'Helvetica,12'
    plot sprintf('%s/trajectories_%d.dat', TABLES, row) using 1:2 with lines lw 0.5 lc rgb pale[row], \
         0 with lines dt 2 lw 1.7 lc rgb 'black'
    unset label 1
    if (row == 1) {
        unset arrow 100
        unset arrow 101
        unset label 100
        unset label 101
    }
    do for [column=1:4] {
        file = sprintf('%s/density_%d_%d.dat', TABLES, row, column)
        unset logscale y
        set yrange [*:*]
        stats file using 3 nooutput
        mean = STATS_min
        set lmargin at screen left[column]
        set rmargin at screen right[column]
        set xrange [-50:50]
        # Major ticks every 40, minor ticks every 10, on both spines.
        set xtics -40,40,40
        set mxtics 4
        set format x ''
        unset xlabel
        if (row == 4) {
            set format x '%g'
            set xlabel '{/:Italic r} − {/:Italic q}  [h^{-1}Mpc]' offset 0,-0.2
        }
        if (LINEAR_Y) {
            set yrange [0:DENSITY_MAX]
            set ytics autofreq
            set mytics 5
            baseline = 0
        } else {
            set logscale y
            set yrange [0.03:DENSITY_MAX]
            set ytics 10
            set mytics 10
            baseline = 0.03
        }
        set format y ''
        unset ylabel
        if (column == 1) {
            set format y '%g'
            if (LINEAR_Y) {
                set ylabel '{/:Italic p} / {/:Italic p}_{/:Normal ref,max}' offset -1.2,0
            } else {
                set ylabel sprintf('{/:Italic p}_{/:Italic t}({/:Italic r} | {/:Italic q} = %g h^{-1}Mpc)', Q[row]) offset -1.2,0
            }
        }
        unset title
        if (row == 1) {
            set title sprintf('{/:Italic a}({/:Italic t}) = %g', times[column]) offset 0,0.45 font 'Helvetica,12'
        }
        set grid xtics mxtics ytics mytics lc rgb '#d4d4d4' dt 3 lw 0.5
        set arrow 10 from 0,baseline to 0,DENSITY_MAX nohead front dt 2 lw 1 lc rgb '#404040'
        set arrow 11 from mean,baseline to mean,DENSITY_MAX nohead front lw 1 lc rgb colours[row]
        if (BANDS[(row-1)*4+column]) {
            alpha = (row == 4 ? 0.10 : 0.06)
            set style fill transparent solid alpha noborder
            plot sprintf('%s/fill_%d_%d.dat', TABLES, row, column) using 1:2:(baseline) \
                with filledcurves between lc rgb colours[row]
            unset style fill
        }
        # Separate pass ensures the complete outline remains above the fill.
        plot sprintf('%s/model_%d_%d.dat', TABLES, row, column) using 1:2 with lines dt 2 lw 1 lc rgb 'black', \
             file using 1:2 with histeps lw 1 lc rgb colours[row]
        unset arrow 10
        unset arrow 11
        unset grid
    }
}
unset multiplot
set output
