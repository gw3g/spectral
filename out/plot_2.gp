
mn= "00"
load "format.gp"

R(k0,k) = 48*16*pi*(k0-k)*(k0-k)/4.
#R(k0,k) = 12*16*pi

set xl "k0/T"
set yl "192 {/Symbol p} x {/Symbol r} (4k@_-^2)"

set yr [-1.2:.2]
set xr [.005:140]
set log x

set key t l
set grid

set tit "I@_{11020}^{(".mn.")}      {\\@}  k=T"
p diag(2,"1.00","+++",mn)   u 1:($2*R($1,1.0))  w l lt 1 t "(+++)" ,\
  diag(2,"1.00","-++",mn)   u 1:($2*R($1,1.0))  w l lt 2 t "(-++)" ,\
  diag(2,"1.00","+-+",mn)   u 1:($2*R($1,1.0))  w l lt 3 t "(+-+)" ,\
  diag(2,"1.00","--+",mn)   u 1:($2*R($1,1.0))  w l lt 4 t "(--+)"

pause -1

