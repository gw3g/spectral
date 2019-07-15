
s=  "+++"
mn= "00"
load "format.gp"

R(k0,k) = 192*pi*abs(k0*k0-k*k)
#R(k0,k) = 8*96*pi*abs(K2(k0,k))/k0 # 01
#R(k0,k) = 4*96*pi*abs(K2(k0,k))/k0 # 10

set xl "k0/T"
set yl "192 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-6.6:0.6]
set xr [.005:140]
set log x

set key t r
set grid

set tit "I@_{11110}^{(".mn.")}   (".s.")"
p diag(5,"0.10",s,mn)   u 1:($2*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(5,"1.00",s,mn)   u 1:($2*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(5,"10.00",s,mn)  u 1:($2*R($1,10.))  w lp lt 3 t "    10.",\
  diag(5,"0.10",s,mn)   u (cut($1,1)):(($3+$4)*R($1,0.1)) w l ls 1 t "OPE" ,\
  diag(5,"1.00",s,mn)   u (cut($1,2)):(($3+$4)*R($1,1.))  w l ls 2 t " " ,\
  diag(5,"10.00",s,mn)  u (cut($1,11)):(($3+$4)*R($1,10.)) w l ls 3 t " " 

pause -1

