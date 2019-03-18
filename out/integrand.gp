set xl 'x/T'
set yl 'y/T'

set xr [-.1:1.1]
set yr [-.1:1.1]
set dgrid3d 80 80
#set zr [-5:5]

set grid
set key t r

set tit "integrand(x,y)/k_0^2 with  k_0=60, k=10"
sp 'data/test_integrand.dat' u 1:2:3 w lp
pause -1
# reread

