
set xl 'x/T'
set yl 'y/T'

set xr [0.1:10.1]
set yr [-5.5:5.5]

set grid
set key t r

set tit "integrand(x,y)/k_0^2 with  k_0=60, k=10"
p 'data/NLO_rhoV_{k=2.09}.dat' u 1:($2) w lp,\
  'data/NLO_rhoV_{k=1.83}.dat' u 1:($2) w lp,\
  'data/NLO_rhoV_{k=3.67}.dat' u 1:($2) w lp,\
  'data/NLO_rhoV_{k=4.19}.dat' u 1:($2) w lp,\
  'data/NLO_rhoV_{k=5.50}.dat' u 1:($2) w lp,\
  'data/NLO_rhoV_{k=6.28}.dat' u 1:($2) w lp
pause -1
# reread

