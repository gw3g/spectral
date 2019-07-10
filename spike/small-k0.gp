
set xl 'k_0/T'

OPE(k0,k) = 4.*k*k/(27.*(k0*k0-k*k)**2)
f(k0,k) = 3.*(k0*k0-k*k)/( 8.*(4.*pi)**3 )
r(k0,k) = 3*(k0*k0-k*k)/(k*k)
#r(k0,k) = 1.

set xr [0:20]
set yr [-.5:.5]

set grid
set key t r

#set tit "{/Symbol r}_V^{nlo} /(gT)^2 - {/Symbol r}_V^{lo}/(2{/Symbol p}T)^2"
#p 'NLO_rho_{k=0.50}.dat' u 1:(($4+r($1,.5)*($5))) w lp t "k/T=0.5",\
#  'NLO_rho_{k=1.00}.dat' u 1:(($4+r($1,1.)*($5))) w lp t "   =1.0",\
#  'NLO_rho_{k=1.50}.dat' u 1:(($4+r($1,1.5)*($5))) w lp t "   =1.5",\
#  'NLO_rho_{k=0.01}.dat' u 1:(($4+r($1,.005)*($5))) w lp t "   =0",\

set tit "{/Symbol r}_{00}/T^2,  O(g^2)"
p 'NLO_rho_{k=1.83}.3.dat' u 1:(($4+r($1,1.8326)*($5))) w lp t "   =1*kn",\
  'NLO_rho_{k=3.67}.3.dat' u 1:(($4+r($1,3.6652)*($5))) w lp t "   =2*kn",\
  'NLO_rho_{k=5.50}.3.dat' u 1:(($4+r($1,5.4978)*($5))) w lp t "   =3*kn",\
  'NLO_rho_{k=1.83}.4.dat' u 1:(($4+r($1,1.8326)*($5))) w lp t "   =3*kn",\
  OPE(x,1.83) lt 1,\
  OPE(x,3.67) lt 2,\
  OPE(x,5.5) lt 3
pause -1
# reread

