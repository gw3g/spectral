s=  "+--"
load "format.gp"

R(k0,k) = 384*pi

set xl "k0/T"
set yl "384 {/Symbol p} {/Symbol r} / T^2"

set yr [-2.1:1.1]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I@_{11110}^{(*)}      (".s.")"
p diag(7,"0.00",s,"00")   u 1:(($2)*R($1,0.0004))  w lp lt 4 t "k/T=.1"  ,\
  diag(7,"0.10",s,"00")   u 1:(($2)*R($1,0.1))  w lp lt 1 t "k/T=.1"  ,\
  diag(7,"1.00",s,"00")   u 1:(($2)*R($1,1.0))  w lp lt 2 t "    1."  ,\
  diag(7,"10.00",s,"00")  u 1:(($2)*R($1,10.))  w lp lt 3 t "    10.",\
  diag(7,"0.10",s,"00")   u (cut($1,1)):(($3+$4)*R($1,0.1)) w l ls 1 t "OPE" ,\
  diag(7,"1.00",s,"00")   u (cut($1,2)):(($3+$4)*R($1,1.))  w l ls 2 t " " ,\
  diag(7,"10.00",s,"00")  u (cut($1,11)):(($3+$4)*R($1,10.)) w l ls 3 t " " 

pause -1

