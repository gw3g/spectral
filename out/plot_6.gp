r0= 64*pi

set xl 'k0/T'
set yl "64 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-100.1:1.1]
set xr [1.005:140]
set log x

set key t r
set grid

set tit "I_{11111}^{(0,0)}, (+++)"
p  'data/diagram6_+++_00_{k=0.004}.dat'  u 1:($2*r0) w lp lt 1 t "k/T=.0",\
  'data/diagram6_+++_00_{k=0.1}.dat'  u 1:($2*r0) w lp lt 2 t "    .1",\
  'data/diagram6_+++_00_{k=1}.dat'  u 1:($2*r0) w lp lt 3 t "    1.",\
  'data/diagram6_+++_00_{k=10}.dat'  u 1:($2*r0) w lp lt 4 t "    10."


pause -1

