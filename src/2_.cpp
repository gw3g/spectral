#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"

struct rho11020 : master {
  double eval();
  rho11020(int _m, int _n, int _s[3]) : master(_m,_n,_s) {}
};
// function(s) for MAIN
master* _11020(int m, int n, int s[3]) {
  master *R =  new rho11020(m,n,s); return R;
}
master* _10120(int m, int n, int s[3]) {
  int s_new[3]; // swap: s_2 <-> s_3
  s_new[0] = s[0];
  s_new[1] = s[1];
  s_new[2] = s[0]*s[1]*s[2];
  master *R =  new rho11020(m,n,s_new); return R;
}

double rho11020::eval() {
  double res=0.;
  res += f(km,s[1])*f(kp,s[4])/kp;
  res -= f(kp,s[1])*f(km,s[4])/km;
  res *= -K2*.0625*( ((double)s[0])*exp(k0)-1. )/k; // coeff = 1/16
  return res*OOFP*I(n,s[2]);
}


/*--------------------------------------------------------------------*/

struct rho11010 : master {
  struct region {
    master *R;
    double operator ()(double p, double xp) {
      double res;
      res = pow(p/k0,(R->n))*f(p,(R->s)[1])*f(k0-p,(R->s)[4]);
      return (((double)(R->s)[0])*exp(k0)-1.)*res;
    }
    double operator ()(double p) { return (*this)(p,k0-p); };
  };
  double eval();
  region ii;
  map *S; 
  rho11010(int _m, int _n, int _s[3]) : master(_m,_n,_s) {
    ii.R = this;
  }
};
// function for MAIN
master* _11010(int m, int n, int s[3]) {
  master *R =  new rho11010(m,n,s); return R;
}

double rho11010::eval() {
  double res=0.;
    if (k0>k) { 
      S = new Finite<region>(ii,km,kp);
      integrate<map> i(*S);
      res += go(i);
    } else {
      S = new SemiInf<region>(ii,km);
      integrate<map> i1(*S);
      res += go(i1);
      S = new SemiInf<region>(ii,kp);
      integrate<map> i2(*S);
      res += go(i2);
    }
  return res;
}

