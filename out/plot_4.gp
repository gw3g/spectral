r0= 64*pi

set xl 'k0/T'
set yl "64 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-0.1:2.1]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I_{11100}^{(0,0)}, (+++)"
p 'data/diagram4_+++_00_{k=0.1}.dat'  u 1:($2*r0) w lp lt 1 t "    .1",\
  'data/diagram4_+++_00_{k=1}.dat'  u 1:($2*r0) w lp lt 2 t "    1.",\
  'data/diagram4_+++_00_{k=10}.dat'  u 1:($2*r0) w lp lt 3 t "    10."

pause -1

