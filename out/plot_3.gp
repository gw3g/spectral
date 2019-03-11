r0= 1. #16*pi*pi

# to remove the LC divergence, multiply by..
fac(k0,k) = 1. #abs((k0-k)/(k0+k))

set xl 'k0/T'
set yl "192 {/Symbol p} x {/Symbol r} k_-/k_+"

set yr [-.05:.1]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I_{11011}^{(0,0)}, (s_0s_1s_2)"
p 'data/diagram3_+++_{k=1}.dat'  u 1:(fac($1,1)*$2/r0) w l lt 1 t "(+++)"

pause -1

