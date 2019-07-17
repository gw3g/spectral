
s=  "+--"
mn= "20"
load "format.gp"

#R(k0,k) = 32*pi*K2(k0,k)*abs(K2(k0,k)) # 00
#R(k0,k) = 32*pi*K2(k0,k)*abs(K2(k0,k))/k0 # 10
#R(k0,k) = 32*pi*(k0*k0-k*k)**2/k0
#R(k0,k) = 32*pi*(k0*k0-k*k) # 11
#R(k0,k) = 32*pi*(k0*k0-k*k)*12./5. # 20, (+++)
R(k0,k) = 32*pi*(k0*k0-k*k)*24 # 20, (+--)
#R(k0,k) = 32*pi
#R(k0,k) = 1.

set xl "k0/T"
set yl "(12/5) x 32 {/Symbol p} K^2 x {/Symbol r} / T^2"

set yr [-2.6:.6]
set xr [.005:140]
set log x

set key b r
set grid

set tit "I@_{11111}^{(".mn.")}   (".s.")"
p diag(6,"0.00",s,mn)   u 1:($2*R($1,0.0004))  w lp lt 4 t "k/T=.0"  ,\
  diag(6,"0.10",s,mn)   u 1:($2*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(6,"1.00",s,mn)   u 1:($2*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(6,"10.00",s,mn)  u 1:($2*R($1,10.))  w lp lt 3 t "    10.",\
  diag(6,"0.10",s,mn)   u (cut($1,1)):(($3+$4)*R($1,0.1)) w l ls 1 t "OPE" ,\
  diag(6,"1.00",s,mn)   u (cut($1,2)):(($3+$4)*R($1,1.))  w l ls 2 t " " ,\
  diag(6,"10.00",s,mn)  u (cut($1,11)):(($3+$4)*R($1,10.)) w l ls 3 t " " 

pause -1

