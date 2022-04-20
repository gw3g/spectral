
s=  "+--"
mn= "00"
mot= "0.00"
load "format.gp"

R(k0,k) = 48*pi*(k0*k0-k*k)
#R(k0,k) = 16*48*pi*(k0*k0-k*k)/k0

set xl "k0/T"
set yl "48 {/Symbol p} K^2 x {/Symbol r}/T^2"

set yr [-.1:1.6]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I@_{11011}^{(".mn.")}   (".s.")"
p diag(3,"0.10",s,mn,mot)   u 1:($2*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(3,"1.00",s,mn,mot)   u 1:($2*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(3,"10.00",s,mn,mot)  u 1:($2*R($1,10.))  w lp lt 3 t "    10." ,\
  diag(3,"0.10",s,mn,mot)   u (cut($1,1)):(($3+$4)*R($1,0.1)) w l ls 1 t "OPE" ,\
  diag(3,"1.00",s,mn,mot)   u (cut($1,2)):(($3+$4)*R($1,1.))  w l ls 2 t " " ,\
  diag(3,"10.00",s,mn,mot)  u (cut($1,11)):(($3+$4)*R($1,10.)) w l ls 3 t " " 

pause -1

