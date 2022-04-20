
s=  "+++"
mn= "00"
mot= "0.00"
load "format.gp"

R(k0,k) = 64*pi

set xl "k0/T"
set yl "64 {/Symbol p} x {/Symbol r}/T^2"

set yr [-0.1:2.1]
set xr [.005:140]
set log x

set key t r
set grid

set tit "I@_{11100}^{(".mn.")}   (".s.")"
p diag(4,"0.10",s,mn,mot)   u 1:($2*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(4,"1.00",s,mn,mot)   u 1:($2*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(4,"10.00",s,mn,mot)  u 1:($2*R($1,10.))  w lp lt 3 t "    10.",\
  diag(4,"0.10",s,mn,mot)   u (cut($1,1)):(($3+$4)*R($1,0.1)) w l ls 1 t "OPE" ,\
  diag(4,"1.00",s,mn,mot)   u (cut($1,2)):(($3+$4)*R($1,1.))  w l ls 2 t " " ,\
  diag(4,"10.00",s,mn,mot)  u (cut($1,11)):(($3+$4)*R($1,10.)) w l ls 3 t " " 

pause -1

