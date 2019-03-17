r0= 2./(pi*pi)

# to remove the LC divergence, multiply by..
fac(k0,k) = 1. #((k0-k)*(k0+k))

set xl 'k0/T'
set yl "48 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-10.1:1.6]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I_{11011}^{(0,0)}, (+++)"
p 'data/diagram4_+++_00_{k=0.1}.dat'  u 1:($2*r0) w lp lt 1 t "    .1",\
  'data/diagram4_+++_00_{k=1}.dat'  u 1:($2*r0) w lp lt 2 t "    1.",\
  'data/diagram4_+++_00_{k=10}.dat'  u 1:($2*r0) w lp lt 3 t "    10."

pause -1

