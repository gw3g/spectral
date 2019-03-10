
#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"

template<typename F>
double in(F *func) {
  map *S; 
  double res=0.;
    if (k0>k) { 
      S = new Finite<F>(*func,km,kp);
      integrate<map> i(*S);
      res += go(i);
    } else {
      S = new SemiInf<F>(*func,km);
      integrate<map> i1(*S);
      res += go(i1);
      S = new SemiInf<F>(*func,kp);
      integrate<map> i2(*S);
      res += go(i2);
    }
  return res;
}

struct H0 {
  // integrand for function H^[0]_m
  int nu, sA, sB;
  double operator ()(double p, double ) {
    double res;
    res = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
    return (((double)sA*sB)*exp(k0)-1.)*res/k;
  };
  double operator ()(double p) { return (*this)(p,k0-p); };
  H0(int _nu, int _sA, int _sB)
    : nu(_nu), sA(_sA), sB(_sB) {}
};


struct H1 {
  // integrand for function H^[1]_m
  // O(eps) terms
  int nu, sA, sB;
  double operator ()(double p, double xp) {
    double res;
    res = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
    if (p<.5*k0) { res *= lga(k*k/( (p-kp)*xp )); }
    else      { res *= lga(k*k/( (p-km)*xp )); }
    return (((double)sA*sB)*exp(k0)-1.)*res/k;
  };
  double operator ()(double p) { return (*this)(p,k0-p); };
  H1(int _nu, int _sA, int _sB)
    : nu(_nu), sA(_sA), sB(_sB) {}
};

struct G0 {
  // integrand for G^[0]_m
  int nu, sA, sB;
  double operator ()(double p, double xp) {
    double res;
    res = f(p,sA)+f(k0-p,sB);
    if (fabs(p-kp)<.1*k) { res *= lga( (p-km)*(p+kp)/(xp*(km+p)) ); }
    else if (( (fabs(p-km)<.1*k)&&(p>fabs(km)/2.) ))
        { res *= lga( xp*(p+kp)/((p-kp)*(km+p)) ); }
    else {return (*this)(p); }
    return (((double)sA*sB)*exp(k0)-1.)*res/k;
  };
  double operator ()(double p) {
    double res;
    res = f(p,sA)+f(k0-p,sB);
    res *= lga( (p-km)*(p+kp)/((p-kp)*(km+p)) ); 
    return (((double)sA*sB)*exp(k0)-1.)*res/k;
  };
  G0(int _nu, int _sA, int _sB)
    : nu(_nu), sA(_sA), sB(_sB) {}
};

/*--------------------------------------------------------------------*/

struct rho11011 : master {
  double eval();
  H0 *h0m, *h0n;
  H1 *h1m, *h1n;
  G0 *g0m;
  map *S; 
  rho11011(int _m, int _n, int _s[3]) : master(_m,_n,_s) {
    h0m = new H0(m,s[1],s[4]);
    h0n = new H0(n,s[2],s[5]);
    h1m = new H1(m,s[1],s[4]);
    h1n = new H1(n,s[2],s[5]);
    g0m = new G0(m,s[1],s[4]);
  }
};
// function for MAIN
master* _11011(int m, int n, int s[3]) {
  master *R =  new rho11011(m,n,s); return R;
}

double rho11011::eval() {
  double res=0.;
  //res += in(h1m);
    if (k0>k) { 
      S = new Finite<G0>(*g0m,0,km);
      integrate<map> i1(*S);
      res += go(i1);
      S = new Finite<G0>(*g0m,km,kp);
      integrate<map> i2(*S);
      res += go(i2);
      S = new SemiInf<G0>(*g0m,kp);
      integrate<map> i3(*S);
      res += go(i3);
    } else {
      S = new Finite<G0>(*g0m,0,kp);
      integrate<map> i1(*S);
      res += go(i1);
      S = new SemiInf<G0>(*g0m,kp);
      integrate<map> i2(*S);
      res += go(i2);
    }
  return res;
  return res;
}
