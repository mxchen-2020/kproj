set term postscript enhanced color 
#set term png truecolor transparent 
set output "plotbnds_up.eps"
#set output "plotbnds_up.png"
set multiplot
set pm3d map
set palette rgb 21,22,23;
set lmargin 13
set rmargin 0.3
#unset ytics
set tmargin 0
set bmargin 0
#set tics out
unset xtics
unset xlabel
set parametric
const=0
set trange [ -21.00000 :   5.00000]
set xrange [0 :     3.177392]
set yrange [ -21.00000 :   5.00000]
set cbrange [0:  0.36]
 
# tics pointing out
set tics out
 
# hide tics on x
unset xtics
set ytics font "Helvetica,32" nomirror # hide tics on the opposite side
#set ytics  -21.000   1.000   5.000 font "Helvetica,32" nomirror # hide tics on the opposite side
 
# set tics width
set border lw 2 
set ylabel offset -5.5,0 
set ylabel "E (eV)" font "Helvetica,32" 
set tmargin 1
set bmargin 1
 
#set label "{/Symbol G}"   font "Helvetica,32" at 0.50,-15.3
set label "K" font "Helvetica,32" at  -0.07943, -22.30000
set label "K" font "Helvetica,32" at   1.39520, -22.30000
set label "K" font "Helvetica,32" at   3.09796, -22.30000
 
splot  "band_kproj_1.dat" 
 
