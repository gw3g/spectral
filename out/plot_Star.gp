r0= 384*pi

set xl 'k0/T'
set yl "64 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-2.1:1.1]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I_{11110}^{(*)}, (+++)"
p 'data/diag.7{k=0.004}.(+++).00.dat'   u 1:($2*r0) w lp lt 1 t "k/T=.0",\
  'data/diag.7{k=0.1}.(+--).00.dat'     u 1:($2*r0) w lp lt 2 t "    .1",\
  'data/diag.7{k=1}.(+--).00.dat'       u 1:($2*r0) w lp lt 3 t "    1.",\
  'data/diag.7{k=10}.(+--).00.dat'      u 1:($2*r0) w lp lt 4 t "    10."


pause -1

