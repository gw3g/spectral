X="00"
#r(k0,k)=k0*k0-k*k
r(k0,k)=1

set xl 'k_0/T'

set xr [0.00:2]
set yr [-.2:.8]

set grid
set key t l

set tit "(K^2/k^2)[ {/Symbol r}@_{".X."}^{nlo} /(gT)^2 - {/Symbol r}@_{".X."}^{lo}  /(2{/Symbol p}T)^2 ]"
p "NLO_rho".X."_{k=0.00}.dat" u 1:(r($1,0)*($3-$2/(2*pi)**2)) w lp t "k/T=0.3",\
  "NLO_rho".X."_{k=0.50}.dat" u 1:(r($1,0.5)*($3-$2/(2*pi)**2)) w lp t "   =1.5",\
  "NLO_rho".X."_{k=1.00}.dat" u 1:(r($1,1)*($3-$2/(2*pi)**2)) w lp t "   =3.0",\
  "NLO_rho".X."_{k=1.50}.dat" u 1:(r($1,1.5)*($3-$2/(2*pi)**2)) w lp t "   =6.0"
pause -1
# reread

