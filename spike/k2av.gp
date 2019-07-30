# cf. Fig.5 in 1310.0164

set yr [-.005:.005]
set xr [.07:150]
set xl "M/T"
set log x

set grid
set key t l

n = -1/(24*16*pi)

set tit "{/Symbol r}/T^2,   w/  k@^2_{av}(M) "
p (4*n) lt 2 not,\
  (n) lt 3 not,\
  (-n) lt 4 not,\
  (n/2) lt 7 not,\
  (-2*n) lt 5 not,\
  (2*n) lt 6 not,\
  'list_k2av.dat'   u 1:2 pt 2 t  "11010 (b:---)",\
  'list_k2av.dat'   u 1:3 pt 4 t  "10110 (b:--+)",\
  'list_k2av.dat'   u 1:4 pt 5 t  "11020 (d:---)",\
  'list_k2av.dat'   u 1:5 pt 6 t  "10120 (d:--+)",\
  'list_k2av.dat'   u 1:6 pt 1 t  "11011 (g:--+)",\
  'list_k2av.dat'   u 1:7 pt 7 t  "{/ZapfDingbats H}  (h':--+)",\
  'list_k2av.dat'   u 1:8 pt 8 t  "11110 (h:--+)",\
  'list_k2av.dat'   u 1:9 pt 9 t  "11110 (m=0,n=1)",\
  'list_k2av.dat'   u 1:10 pt 10 t  "11111 (j:--+)",\
  'list_k2av.dat'   u 1:11 pt 11 t  "11111 (m=2)"

pause -1
