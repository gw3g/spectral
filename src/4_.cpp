
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
struct region {
  Master *R;
  int _s[3]={(R->s)[0],(R->s)[1],(R->s)[2]};
  int _m = R->m;
  int _n = R->n;
  virtual double func(double,double) = 0;
  double operator ()(double x,double y) {return func(x,y);};
};

template<typename F>
Finite<F> make(F f,double a, double b) {
  return Finite<F>(f,a,b);
}

template<typename F>
SemiInf<F> make(F f,double a) {
  return SemiInf<F>(f,a);
}
/*
struct above_LC_1A : region {
  double func(double x, double y) {
  return make([&](double p, double pd) { 
  return make([&](double q, double qd) {
      double pqr[3] = {p,q,k0-p-q};
      return W_iv(p,q)*fff(pqr,_m,_n,_s);
  },//  q_min,  q_max
        kp-p,   p     )(y);
  },//  p_min,  p_max
        .5*kp,  km    )(x);
  };
};*/

struct rho11100 : Master {
  double eval();
  double above_LC_1A(double,double);
  //double above_LC_2A(double,double);

  rho11100(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
    //reg = new above_LC_1A();
    //(*reg).R = this;
    //above_LC_1A.R = this;
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
  //res += integrate_2d(i);
  //prt = above_LC_1A;
  //res += integrate_2d(above_LC_1A);
  //res += above_LC_1A(.4,.3);
  //res += above_LC_2A(.2,.04);
  //res += above_LC_1A(.4,.3);
  return res;
}


double rho11100::above_LC_1A(double x, double y) 
{
  int _s[3]={(this->s)[0],(this->s)[1],(this->s)[2]};
  int _m = this->m;
  int _n = this->n;
  return make([&](double p, double pd) { 
  return make([&](double q, double qd) {
      double pqr[3] = {p,q,k0-p-q};
      return W_iv(p,q)*fff(pqr,_m,_n,_s);
  },//  q_min,  q_max
        kp-p,   p     )(y);
  },//  p_min,  p_max
        .5*kp,  km    )(x);
}
