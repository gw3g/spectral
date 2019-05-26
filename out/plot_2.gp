
s=  "+++"
mn= "00"
load "format.gp"

R(k0,k) = 12*16*pi*(k0-k)*(k0-k)/4.

set xl "k0/T"
set yl "192 {/Symbol p} x {/Symbol r} (4k_-^2)"

set yr [-3.6:.6]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I_{11020}^{(".mn.")}, (".s.")"
p diag(2,"0.10",s,mn)   u 1:($2*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(2,"1.00",s,mn)   u 1:($2*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(2,"10.00",s,mn)  u 1:($2*R($1,10.))  w lp lt 3 t "    10."

pause -1

