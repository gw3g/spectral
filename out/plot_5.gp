r0= 192*pi

set xl 'k0/T'
set yl "192 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-0.1:2.1]
set xr [.005:140]
set log x

set key t r
set grid

set tit "I_{11110}^{(0,0)}, (+++)"
p 'data/diagram4_+++_00_{k=0.004}.dat'  u 1:($2*r0) w lp lt 4 t "k/T=.0",\
  'data/diag.5{k=0.1}.(+++).00.dat'  u 1:($2*r0) w lp lt 1 t "    .1",\
  'data/diag.5{k=1}.(+++).00.dat'  u 1:($2*r0) w lp lt 2 t "    1.",\
  'data/diag.5{k=10}.(+++).00.dat'  u 1:($2*r0) w lp lt 3 t "    10."

pause -1

