#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"

double F(int nu, int sA, int sB) {
  // here F_m for m={0,1,2}
  double res=0.;
  res += f(km,sA)*f(kp,sB)/kp;
  res -= f(kp,sA)*f(km,sB)/km;
  res *= ((double)sA*sB)*exp(k0)-1.; // = F_0
  if (nu>0) {
    res *= k0;
    res += ( 2.-exp(kp)*((double)(sA+sB)) )*f(kp,sA)*f(kp,sB);
    res += sgn(km)*( f(fabs(km),sA)+f(fabs(km),sB) ); // = F_1
  }
  if (nu>1) {
    res *= k0;
    res -= (km>0)?k:0.;
    res += 2.*(( kp-fabs(km)-lga( (exp(kp)-sB)/(exp(fabs(km))-sB) ) ));
    res += f(fabs(km),sB)*( k0 + f(fabs(km),sA)*exp(fabs(km))*km*((double)(sB-sA)) )*sgn(km);
    res -= f(fabs(kp),sB)*( k0 + f(fabs(kp),sA)*exp(fabs(kp))*kp*((double)(sB-sA)) ); // = F_2
  }

  res *= -.0625*OOFP/k; // coeff = 1/16
  return res;
}

/*--------------------------------------------------------------------*/

struct rho11020 : Master {
  double eval();
  rho11020(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=2; }
};
// function(s) for MAIN
Master* _11020(int m, int n, int s[3]) {
  Master *R =  new rho11020(m,n,s); return R;
}
Master* _10120(int m, int n, int s[3]) {
  int s_new[3]; // swap: s_2 <-> s_3
  s_new[0] = s[0];
  s_new[1] = s[1];
  s_new[2] = s[0]*s[1]*s[2];
  Master *R =  new rho11020(m,n,s_new); return R;
}

double rho11020::eval() { return F(m,s[1],s[4])*I(n,s[2]); }


/*--------------------------------------------------------------------*/

struct rho11010 : Master {
  struct region {
    Master *R;
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
  rho11010(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
    ii.R = this;
  }
};
// function for MAIN
Master* _11010(int m, int n, int s[3]) {
  Master *R =  new rho11010(m,n,s); return R;
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

