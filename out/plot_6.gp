r0= 32*pi #*12./5.
#r0 = 32*(4*pi)**3
#r0=1.

set xl 'k0/T'
set yl "32 {/Symbol p} K^2 x {/Symbol r} / T^2"
#set yl "( {/Symbol r}-OPE(T^2) )/OPE(T^4)"

set yr [-2.1:2.1]
set xr [.005:140]
set log x

set key b r
set grid

set tit "I_{11111}^{(1,1)}, (+++)"
p 'data/diag.6{k=0.004}.(+++).00.dat'   u 1:($2*r0) w lp lt 1 t "k/T=.0",\
  'data/diag.6{k=0.10}.(+--).11.dat'     u 1:($2*r0) w lp lt 2 t "    .1",\
  'data/diag.6{k=1.00}.(+--).11.dat'       u 1:($2*r0) w lp lt 3 t "    1.",\
  'data/diag.6{k=10.00}.(+--).11.dat'      u 1:($2*r0) w lp lt 4 t "    10."
pause -1

