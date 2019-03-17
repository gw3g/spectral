
#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"

//W_IV
double W_iv(double p, double q) {
  double r=k0-p-q;
  return (  fabs(p-kp)-fabs(p-km)+
            fabs(q-kp)-fabs(q-km)+
            fabs(r-kp)-fabs(r-km)-min(k0,k) )/(2.*k);
}

/*--------------------------------------------------------------------*/

struct rho11100 : Master {
  double eval();
  double integrand(double,double); // supported on [0,1]x[0,1]

  struct inner {
    rho11100 *R;
    double _x; // x-dependence "stands by"
    double operator ()(double y) { return (R->integrand)(_x,y); }
  };
  struct outer {
    inner f2;
    double operator ()(double x) {
      f2._x = x;
      integrate<inner> I(f2); // do the y-integral
      return go(I);
    };
  };
  double quad() {
    outer f1;
    f1.f2.R = this;
    integrate<outer> I(f1); // do the x-integral
    return go(I)-K2/8.;
  }

  rho11100(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
  }
};
// function for MAIN
Master* _11100(int m, int n, int s[3]) {
  Master *R =  new rho11100(m,n,s); return R;
}

/*--------------------------------------------------------------------*/

double rho11100::eval() {
  double res=0.;
  //double p = 1.9, q= 1.3;
  //double pqr[3] = {p,q,k0-p-q};
  //int cut[3] = {s[0],s[1],s[2]};
  //res += fff(pqr,m,n,cut);
  //return -res*.25*pow(OOFP,3);
  //res += integrate_2d(i);
  //prt = above_LC_1A;
  res += quad();
  //res += quad((*reg).func);
  //res += (*reg).func(.4,.3);
  //res += (double)reg.mm;
  //res += above_LC_2A(.2,.04);
  //res += above_LC_1A(.4,.3);
  //res += (*reg)();
  return res;
}


double rho11100::integrand(double x, double y) 
{
  int _s[3]={(this->s)[0],(this->s)[1],(this->s)[2]};
  int _m = this->m;
  int _n = this->n;
  double res=0.; 
  /*res+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double r = k0-p-q; r= (fabs(r)<min(k,fabs(km))*1e-2) ? qd : r;
      double pqr[3] = {p,q,r};
      return W_iv(p,q)*fff(pqr,_m,_n,_s);
  },  max(0.,km-p), km  )(y); },  0.,kp    )(x);
  res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double r = k0-p-q; r= (fabs(r)<min(k,fabs(km))*1e-2) ? qd : r;
      double pqr[3] = {p,q,r};
      return W_iv(p,q)*fff(pqr,_m,_n,_s);
  },  km-q,   kp     )(y); }, 0.,  km    )(x);*/
  /*res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      q= (fabs(q)<k*1e-2) ? pd : q;
      double pqr[3] = {p,q,k0-p-q};
      return 3.*W_iv(p,q)*fff(pqr,_m,_n,_s);
  },  km, kp-q     )(y); }, -k0, .0    )(x);*/
  res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double pqr[3] = {p,q,k0-p-q};
      return 3.*W_iv(p,q)*fff(pqr,_m,_n,_s);
  },  km, kp-q     )(y); }, -km,km    )(x);
  return res;
}
