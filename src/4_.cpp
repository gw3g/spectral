
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
  static double above_LC_1A(double,double);
  rho11100(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
  }
};
// function for MAIN
Master* _11100(int m, int n, int s[3]) {
  Master *R =  new rho11100(m,n,s); return R;
}

double rho11100::eval() {
  double res=0.;
  //double p = 1.9, q= 1.3;
  //double pqr[3] = {p,q,k0-p-q};
  //int cut[3] = {s[0],s[1],s[2]};
  //res += fff(pqr,m,n,cut);
  //return -res*.25*pow(OOFP,3);
  res += integrate_2d(above_LC_1A);
  //res += above_LC_1A(.4,.3);
  return res;
}

template<typename F>
Finite<F> make(F f,double a, double b) {
  return Finite<F>(f,a,b);
}

double rho11100::above_LC_1A(double x, double y) 
{
  //int _s[3]={(this->s)[0],(this->s)[1],(this->s)[2]};
  int _s[3]={1,1,1};
  //int _m = this->m;
  //int _n = this->n;
  return make([&](double p, double pd) { 
  return make([&](double q, double qd) {
      double pqr[3] = {p,q,k0-p-q};
      return W_iv(p,q)*fff(pqr,0,0,_s);
  },//  q_min,  q_max
        kp-p,   p     )(y);
  },//  p_min,  p_max
        .5*kp,  km    )(x);
}
