r0=12*16*pi

# to remove the LC divergence, multiply by..
fac(k0,k) = ((k0-k)/(k0+k))

set xl 'k0/T'
set yl "192 {/Symbol p} x {/Symbol r} k_-/k_+"

set yr [-.1:1.1]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I_{11020}^{(0,0)}, (s_0s_1s_2)"
p 'data/diagram2_+++_{k=1}.dat'  u 1:($2*r0*fac($1,1)) w l lt 1 t "(+++)",\
  'data/diagram2_-++_{k=1}.dat'  u 1:($2*r0*fac($1,1)) w l lt 2 t "(-++)",\
  'data/diagram2_+-+_{k=1}.dat'  u 1:($2*r0*fac($1,1)) w l lt 3 t "(+-+)",\
  'data/diagram2_--+_{k=1}.dat'  u 1:($2*r0*fac($1,1)) w l lt 4 t "(--+)"

pause -1

