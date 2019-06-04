
set xl 'k_0/T'

f(k0,k) = 3.*(k0*k0-k*k)/( 8.*(4.*pi)**3 )
r(k0,k) = 1. #(k0*k0-k*k)/(k*k)
ope(k0,k)= 1*(4*pi/27)*k*k/(k0*k0-k*k)

set log x
set xr [0.01:100]
set yr [-1:1]

set grid
set key t r

set tit "{/Symbol r}_V^{nlo} /(gT)^2 - {/Symbol r}_V^{lo}/(2{/Symbol p}T)^2"
p 'NLO_rhoV_{k=0.30}.dat' u 1:(r($1,.3)*($3-$2/(2*pi)**2)) w lp t "k/T=0.3",\
  'NLO_rhoV_{k=1.50}.dat' u 1:(r($1,1.5)*($3-$2/(2*pi)**2)) w lp t "   =1.5",\
  'NLO_rhoV_{k=3.00}.dat' u 1:(r($1,3.)*($3-$2/(2*pi)**2)) w lp t "   =3.0",\
  'NLO_rhoV_{k=6.00}.dat' u 1:(r($1,6.)*($3-$2/(2*pi)**2)) w lp t "   =6.0",\
  'NLO_rhoV_{k=9.00}.dat' u 1:(r($1,9.)*($3-$2/(2*pi)**2)) w lp t "   =9.0",\
  ope(x,.3) lt 1 not,\
  ope(x,1.5) lt 2 not
pause -1
# reread

