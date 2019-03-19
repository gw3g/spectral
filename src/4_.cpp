
#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"

//W_IV
double W_iv(double p, double q) {
  double r=k0-p-q;
  return (  fabs(p-kp)-fabs(p-km)+
            fabs(q-kp)-fabs(q-km)+
            fabs(r-kp)-fabs(r-km)-min(k0,k) )/(2.*k); }

/*--------------------------------------------------------------------*/

struct rho11100 : Master {
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
  double eval() {
    double res;
    outer f1;
    f1.f2.R = this;
    integrate<outer> I(f1); // do the x-integral
    res += go(I) + ( (k0>k) ? -K2/8. : 0. );
    return res*pow(OOFP,3);
  }
rho11100(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
  }
};
// function for MAIN
Master* _11100(int m, int n, int s[3]) {
  Master *R =  new rho11100(m,n,s); return R;
}

#include <fstream>
using namespace std;
void print_integrand(int m, int n, int s[3]) {
  ofstream fout;
  rho11100 *R =  new rho11100(m,n,s);
  fout.open("out/data/test_integrand.dat");
  double x=-.0125*.5;
  double y;
  for (int i=0;i<80;i++) {
    x+=.025*.5;
    y=-.0125*.5;
    for (int j=0;j<80;j++) {
      y+=.025*.5;
      fout << x << "    " << y << "    " << (  (R->integrand)(x,y)  ) << endl; }
  }
  fout.close(); 
}

/*--------------------------------------------------------------------*/

double rho11100::integrand(double x, double y) 
{
  int _s[3]={(this->s)[0],(this->s)[1],(this->s)[2]};
  int _m = this->m;
  int _n = this->n;
  double res=0.; 
  //double pqr[3];
  if (k0>k) {

  res+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      return W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
  },  0., k0-p  )(y); },  km,kp    )(x);

  res+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      return W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
  },  km-p,   kp     )(y); }, 0.,  km    )(x);

  res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  km, kp-q     )(y); }, -k0, .0    )(x);

  res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  km, kp-q     )(y); }, -k0    )(x);

  res+= // pq := p-q,
   make([&](double r, double rd) { return make([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*fff(p,q,r,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,r,_m,_n,_s);
      return .5*temp;
  },  0.,k0-2.*km-r )(y); }, -k0,0.    )(x); 
  res+= // pq := p-q,
   make([&](double r, double rd) { return make([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*fff(p,q,r,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,r,_m,_n,_s);
      return .5*temp;
  },  0.,k0-2.*km-r )(y); }, -k0    )(x); 

  } else { // Below the LC:

  res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  km, kp-q     )(y); }, km    )(x);

  res+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0.;
      temp += W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  kp-p, 0.     )(y); }, kp    )(x);

  res+= // pq := p-q,
   make([&](double r, double rd) { return make([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*fff(p,q,r,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,r,_m,_n,_s);
      return .5*temp;
  },  0.,kp )(y); }, km    )(x); 

  res+= // pq := p-q,
   make([&](double r, double rd) { return make([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*fff(p,q,r,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,r,_m,_n,_s);
      return .5*temp;
  },  km, 0. )(y); }, kp    )(x); 

  res+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0.;
      temp += W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  0.,p-kp     )(y); }, kp    )(x);

  res+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_iv(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  q-km, 0.     )(y); }, km    )(x);



  }
  return res;
}
