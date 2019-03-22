

#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"

//W_VI
double W_vi(double p, double q) {
  double r=k0-p-q,
         c_p=clp(p,km,kp),
         c_q=clp(q,km,kp),
         c_r=clp(r,km,kp);

  // kinematic exclusions:
  if (k0>k) {
              if ( (p>kp)&&(q<km)&&(q>kp-p) ) return 0.;
              if ( (q>kp)&&(p<km)&&(p>kp-q) ) return 0.;
              if ( (p<km)&&(q<km)&&(q+p<km) ) return 0.;
  }

  double res; bool reg_p = (km>p)||(kp<p), reg_q = (km>q)||(kp<q); // reg_p=0 if p in [km,kp]

  double L1= lga( (k0-p-c_r)*(k0-q-c_r)/( (c_q-p)*(c_p-q) ) ),
         L2= reg_q ? lga( (q-km)/(q-c_q) ) : (  lga( ((p-q)*(km-q))/(p*q) ) ),
         L3= reg_p ? lga( (p-km)/(p-c_p) ) : (  lga( ((p-q)*(km-p))/(p*q) ) ),
         L6= lga( (k0-p-c_r)*(k0-q-c_r)/( p*q ) ),
         L7= reg_p ? lga( p*q*(k0-q-c_p)/( (km-q)*(kp-q)*(p-c_p) ) ) : 0.,
         L8= reg_q ? lga( p*q*(k0-p-c_q)/( (km-p)*(kp-p)*(q-c_q) )  ) : 0.;

  if (k0>k) { res =L1; res+=L2+L3; }
  else      { res =L6; res+=L7+L8; }

  return res;
}

/*--------------------------------------------------------------------*/

struct rho11111 : Master {
  double integrand(double,double); // supported on [0,1]x[0,1]

  struct inner {
    rho11111 *R;
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
    res += go(I);
    return res*pow(OOFP,3);
  }

  rho11111(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
  }
};
// function for MAIN
Master* _11111(int m, int n, int s[3]) {
  Master *R =  new rho11111(m,n,s); return R;
}

/*--------------------------------------------------------------------*/

double rho11111::integrand(double x, double y) { 
  int _s[3]={(this->s)[0],(this->s)[1],(this->s)[2]};
  int _m = this->m;
  int _n = this->n;
  double res=0.; 

  if (k0>k) {

  res += // 1A
   make([&](double p, double pd) { return make([&](double q, double qd) {
         double temp = 0.;
         if (k0>3.*k) {
          temp +=  W_vi(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
          temp +=  W_vi(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
          } 
          return temp;
  },  kp-p, p  )(y); },  .5*kp, km    )(x);//*/

  res += // 1B
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0.;
      double pm = km-p, qm = km-q;
      if (pm<1e-1*( min(km,kp-q)-max(km-q,q) )) { pm = pd; }
      temp += lga(p*q/(kp*qm))*fff(p,q,k0-p-q,_m,_n,_s)/(k0-p-q);
      temp += lga(p*q/(kp*qm))*fff(q,p,k0-p-q,_m,_n,_s)/(k0-p-q);
      return temp;
  },  max(km-q,q), min(km,kp-q) )(y); }, 0.,min(km,.5*kp)   )(x);//*/

  res+= // 1D, TODO
   make([&](double p, double pd) { return make([&](double q, double qd) {
         return 0.;
  },  kp-p, km  )(y); },  max(k,km), kp    )(x);

  res+= // 2A
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      temp += W_vi(p,q)*fff(p,q,r,_m,_n,_s)/r;
      temp += W_vi(q,p)*fff(q,p,r,_m,_n,_s)/r;
      return temp;
  },  km-p  )(y); },  kp )(x); //*/

  res+= // 2B
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      double pp = kp-p, qp = kp-q;
      if (qp<1e-1*k) { qp = pd; }
      //temp += W_vi(p,q)*fff(p,q,r,_m,_n,_s)/r;
      //temp += W_vi(q,p)*fff(q,p,r,_m,_n,_s)/r;
      temp += lga(p*q/(pp*qp))*fff(p,q,r,_m,_n,_s)/r;
      temp += lga(p*q/(pp*qp))*fff(q,p,r,_m,_n,_s)/r;
      return temp;
  },  km-p,kp-p  )(y); },  kp )(x); //*


  /*res+= // 2C
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      temp += W_vi(p,q)*fff(p,q,r,_m,_n,_s)/r;
      temp += W_vi(q,p)*fff(q,p,r,_m,_n,_s)/r;
      return temp;
  },  km-p  )(y); },  km,kp )(x); //*/

  res+= // 4A
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double r = k0-p-q;
      return W_vi(p,q)*fff(p,q,r,_m,_n,_s)/r;
  },  kp  )(y); },  kp    )(x); //*/

  /*res+= // 4B?
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      temp += W_vi(p,q)*fff(p,q,r,_m,_n,_s)/r;
      temp += W_vi(q,p)*fff(q,p,r,_m,_n,_s)/r;
      return temp;
  },  km,kp  )(y); },  kp    )(x); //*/

  /*res+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0.;
      temp += W_vi(p,q)*fff(p,q,k0-p-q,_m,_n,_s);
      temp += W_vi(q,p)*fff(q,p,k0-p-q,_m,_n,_s);
      return temp;
  },  kp,km-p  )(y); },  -k    )(x);*/

  } else {
    res += 0.;
  }
  return res;
}
