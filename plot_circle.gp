set term pdf
set output "circle_evol.pdf"

unset key
unset colorbox
unset xtics
unset ytics

set xrange [0:200]
set yrange [0:200]
set multiplot  
set size 0.35, 0.55
set origin 0.1, 0.5
set title "t=0"
plot "Circulo1.dat" matrix with image
set size 0.35, 0.55
set origin 0.5, 0.5
set title "t=400"
plot "Circulo2.dat" matrix with image
set size 0.35, 0.55
set origin 0.1, 0.0
set title "t=800"
plot "Circulo3.dat" matrix with image
set size 0.35, 0.55
set origin 0.5, 0.0
set title "t=1200"
plot "Circulo4.dat" matrix with image
