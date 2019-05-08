
set xl 'k_0/T'

f(k0,k) = 3.*(k0*k0-k*k)/( 8.*(4.*pi)**3 )

set log x
set xr [0.01:100]
set yr [-.1:.1]

set grid
set key t r

set tit "{/Symbol r}_V^{nlo} /(gT)^2 - {/Symbol r}_V^{lo}/(2{/Symbol p}T)^2"
p 'NLO_rhoV_{k=0.30}.dat' u 1:($3-$2/(2*pi)**2) w lp t "k/T=.3",\
  'NLO_rhoV_{k=3.00}.dat' u 1:($3-$2/(2*pi)**2) w lp t "   =3.",\
  'NLO_rhoV_{k=6.00}.dat' u 1:($3-$2/(2*pi)**2) w lp t "   =6.",\
  'NLO_rhoV_{k=9.00}.dat' u 1:($3-$2/(2*pi)**2) w lp t "   =9."
pause -1
# reread

