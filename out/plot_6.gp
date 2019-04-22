r0= 16*pi #*12./5.
#r0 = (4*pi)**3

set xl 'k0/T'
set yl "(12./5.) x 32 {/Symbol p} K^2 x {/Symbol r} / T^2"

set yr [-2.1:2.1]
set xr [.005:140]
set log x

set key b r
set grid

set tit "I_{11111}^{(2,0)}, (+++)"
p 'data/diag.6{k=0.004}.(+++).01.dat'   u 1:($2*r0) w lp lt 1 t "k/T=.0",\
  'data/diag.6{k=0.1}.(+++).11.dat'     u 1:($2*r0) w lp lt 2 t "    .1",\
  'data/diag.6{k=1}.(+++).11.dat'       u 1:($2*r0) w lp lt 3 t "    1.",\
  'data/diag.6{k=10}.(+++).11.dat'      u 1:($2*r0) w lp lt 4 t "    10.",\

pause -1

