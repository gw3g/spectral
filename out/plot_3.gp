
s=  "+++"
mn= "00"
load "format.gp"

R(k0,k) = 48*pi*(k0*k0-k*k)

set xl "k0/T"
set yl "48 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-.1:1.6]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I@_{11011}^{(".mn.")}   (".s.")"
p diag(3,"0.10",s,mn)   u 1:($2*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(3,"1.00",s,mn)   u 1:($2*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(3,"10.00",s,mn)  u 1:($2*R($1,10.))  w lp lt 3 t "    10."

pause -1

