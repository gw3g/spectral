#include "core.hh"
#include "quad.hh"
//#include "trapezoid.hh"
#include "map.hh"

using namespace std;

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
  double F_123(double,double,double);
  double F_345(double,double,double);
  double F_14(double,double);
  double F_25(double,double);

  /*struct inner {
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
  };*/
  double eval() {
    double res, err;
    //outer f1;
    //f1.f2.R = this;
    //integrate<outer> I(f1); // do the x-integral
    //res = go(I);
  double epsabs = 1e-3, epsrel = 0;
size_t limit = 1e5;

  quad wsp1(limit);
  quad wsp2(limit);

  auto outer = make_gsl_function( [&](double x) {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
        return (this->integrand)(x,y);
             //return (this->integrand)(.5*(1.-x),y)+(this->integrand)(1.-.5*x,y);
        } );
    gsl_integration_qag(inner, .0+1e-16,1., epsabs, epsabs, 
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag(  outer, .0+1e-16,1., epsabs, epsabs, 
                         limit, 6, wsp2, &res, &err  );

    return res*CUBE(OOFP);
  }

  rho11111(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
  }
};
// function for MAIN
Master* _11111(int m, int n, int s[3]) {
  Master *R =  new rho11111(m,n,s); return R;
}

#include <fstream>
using namespace std;
void print_integrand(int m, int n, int s[3]) {
  ofstream fout;
  rho11111 *R =  new rho11111(m,n,s);
  fout.open("out/data/test_integrand.dat");
  double x=-.0125*.5;
  double y;
  for (int i=0;i<80;i++) {
    x+=.025*.5;
    y=-.0125*.5;
    for (int j=0;j<80;j++) {
      y+=.025*.5;
      fout << x << "    " << y << "    " << (  
          (R->integrand)(x,y)  
          ) << endl; }
  }
  fout.close(); 
}

/*--------------------------------------------------------------------*/
// all the thermal weights for cuts from this topology

double rho11111::F_123(double p, double q, double r) { 
  double res;

  int s1 = (this->s)[1],
      s2 = (this->s)[2],
      s3 = (this->s)[3];

  int _m = this->m;
  int _n = this->n;

  double fp=f(p,s1), fq=f(q,s2), fr= f(r,s3);

  res = ( ((double) s1*s2*s3)*exp(k0)-1. )*fp*fq*fr;
  res*= pow(k0,-_m-_n)*pow(p,_m)*pow(q,_n) ;

  return res;
}

double rho11111::F_345(double r, double l, double v) {
  double res;

  int s3 = (this->s)[3],
      s4 = (this->s)[4],
      s5 = (this->s)[5];

  int _m = this->m;
  int _n = this->n;

  double fr=f(r,s3), fl=f(l,s4), fv= f(v,s5);

  res = ( ((double) s3*s4*s5)*exp(k0)-1. )*fr*fl*fv;
  res*= pow(k0,-_m-_n)*pow(k0-l,_m)*pow(k0-v,_n) ;

  return res;
}

double rho11111::F_14(double p, double l) {
  double res;

  int s1 = (this->s)[1],
      s4 = (this->s)[4];

  int _m = this->m;

  double fp=f(p,s1), fl=f(l,s4);

  res = ( ((double) s1*s4)*exp(k0)-1. )*fp*fl;
  res*= pow(k0,-_m)*pow(p,_m);

  return res;
}

double rho11111::F_25(double q, double v) {
  double res;

  int s2 = (this->s)[2],
      s5 = (this->s)[5];

  int _n = this->n;

  double fq=f(q,s2), fv=f(v,s5);

  res = ( ((double) s2*s5)*exp(k0)-1. )*fq*fv;
  res*= pow(k0,-_n)*pow(q,_n);

  return res;
}

/*--------------------------------------------------------------------*/

double rho11111::integrand(double x, double y) { 

  // statistics:
  int s0 = (this->s)[0],
      s1 = (this->s)[1],
      s2 = (this->s)[2],
      s3 = (this->s)[3],
      s4 = (this->s)[4],
      s5 = (this->s)[5];

  double res=0.; 

  if (k0>k) {
  double _1A, _1B,      _1D,
         _2A, _2B, _2C, _2D, // 3 by reflection
         _4A, _4B,      _4D;

  //  #[ from '1'

  _1A =
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0., r=k0-p-q;
      double pm = km-p, qm = km-q;
      double pp = kp-p, qp = kp-q;
      if (k0>3.*k) {
        //pm = (fabs(pm)<1e-1*(km-max(kp-q,q))) ? pd : pm;
        temp +=  lga(pp*qp/(pm*qm))*( F_123(p,q,r)+F_123(q,p,r) );
      };
      return temp/r;
  },  max(kp-q,q), km  )(y); },  k, km    )(x);//*/

  _1B =
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, qm = km-q; // p=km (!)
      //pm = (fabs(pm)<1e-1*( min(km,kp-q)-max(km-q,q) )) ? pd : pm;
      temp += lga( p*q/(pm*qm) )*( F_123(p,q,r)+F_123(q,p,r) );
      return temp/r;
  },  max(km-q,q), min(km,kp-q) )(y); }, 0.,min(km,.5*kp)   )(x);//*/

  _1D = 0.; // see 4B

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, qm = km-q;
      double pp = kp-p, qp = kp-q;
      //pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      temp += lga(pm*qm/(pp*qp))*( F_123(p,q,r) + F_123(q,p,r) );
      //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      return temp/r;
  },  km-p  )(y); },  kp )(x); //*/

  _2B =
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      double pp = kp-p, qp = kp-q;
      //pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      temp += lga(p*q/(pp*qp))*( F_123(p,q,r) + F_123(q,p,r) );
      //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      return temp/r;
  },  km-p,kp-p  )(y); },  kp )(x); //*/

  _2C =
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, qm = km-q;
      double pp = kp-p, qp = kp-q;
      pm = (fabs(pm)<1e-1*k) ? pd : pm; // p=km (!)
      pp = (fabs(pp)<1e-1*k) ? pd : pp; // p=kp (!)
      double q_ = (fabs(q)<k0) ? qd : q; // q=0 (!)
      //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp += lga(pm*qm/(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 2.*lga( q*(k0-q)/(qp*qm) )*F_14(p,k0-p)*(.5+f(r,s3)); // virtual
      return temp/r;
  },  km-p  )(y); },  km,kp )(x); //*/

  _2D =
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      double qm = km-q, qp = kp-q;
      double q_ = (fabs(q)<1e-1*fabs(km-p)) ? qd : q; // q=0 (!)
      temp += 2.*lga( q*(k0-q)/(qp*qm) )*F_14(p,k0-p)*(.5+f(r,s3)); // virtual
      return temp/r;
  },  km-p, 0.  )(y); },  km,kp    )(x); //*/

  //  #] from '2','3'\
      #[ from '4'

  _4A =
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double r = k0-p-q;
      return W_vi(p,q)*F_123(p,q,r)/r;
  },  kp  )(y); },  kp    )(x); //*/

  _4B =
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, pp = kp-p;
      temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 2.*lga( p*(k0-p)/(pp*pm) )*F_14(q,k0-q)*(.5+f(r,s3)); // virtual
      return temp/r;
  },  k0  )(y); },  km,kp    )(x);
  if (k0>3.*k) { 
    _4B+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, pp = kp-p;
      double l = k0-p; //l = (fabs(l)<1e-1*(q-km)) ? pd : l;
      temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 2.*lga( p*l/(pp*pm) )*F_14(q,k0-q)*(.5+f(r,s3)); // virtual
      temp += 2.*lga( p*l/(pp*pm) )*F_14(k0-q,q)*(.5+f(r,s3)); // virtual
      return temp/r;
  },  k0+km-q,k0  )(y); },  km,kp    )(x);
    _4B+=
   make([&](double q, double qd) { return make([&](double p, double pd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, pp = kp-p;
      double l = k0-p; //l = (fabs(l)<1e-1*(km-k)) ? pd : l;
      temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp -= W_vi(k0-p,k0-q)*F_123(k0-p,k0-q,-r); // from 1D & 1B
      temp -= W_vi(k0-q,k0-p)*F_123(k0-q,k0-p,-r);
      temp += 2.*lga( p*l/(pp*pm) )*F_14(q,k0-q)*(.5+f(r,s3)); // virtual
      temp += 2.*lga( p*l/(pp*pm) )*F_14(k0-q,q)*(.5+f(r,s3)); // virtual
      return temp/r;
  },  2.*kp-q,k0+km-q  )(y); },  km,kp    )(x)  ;
  } else if (k0<3.*k) { 
    _4B+=
   make([&](double p, double pd) { return make([&](double q, double qd) {
      double temp = 0., r = k0-p-q;
      double pm = km-p, pp = kp-p;
      double l = k0-p; //l = (fabs(l)<1e-1*fabs(k-k0+p)) ? pd : l;
      temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 2.*lga( p*l/(pp*pm) )*F_14(q,k0-q)*(.5+f(r,s3)); // virtual
      temp += 2.*lga( p*l/(pp*pm) )*F_14(k0-q,q)*(.5+f(r,s3)); // virtual
      return temp/r;
  },  k0+km-p, kp  )(y); },  kp,k0    )(x); } //*/

  _4D =
   make([&](double pq, double pqd) { return make([&](double r, double rd) {
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double qm = km-q, qp = kp-q;
      double pm = km-p, pp = kp-p;
      double h = 2.*fabs(k-pq);
      //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
      //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
      //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
      //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;

      temp += F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qp/(q*q) );
      temp += F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pp/(p*p) );
      temp += F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qp/(q*q) );
      temp += F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pp/(p*p) );
      return .5*temp/r;
  },  -k+pq,k-pq  )(y); }, max(0.,k-km),k    )(x); //*/
  _1D += // corner, 1D+B
   make([&](double pq, double pqd) { return make([&](double r, double rd) {
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double qm = km-q, qp = kp-q;
      double pm = km-p, pp = kp-p;
      double h = fabs( min(k,km)-fabs(pq-k) );

      //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
      //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
      //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
      //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

      temp += lga( pp*qp/(q*p) )*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 4.*lga( q*(k0-q)/(qp*qm) )*F_14(p,k0-p)*(.5+f(r,s3)); // virtual
      temp -= lga( pp*qp/((k0-p)*(k0-q)))*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );
      return .5*temp/r;
  },  fabs(pq-k),min(k,km)  )(y); }, k-min(k,km),k+min(k,km)  )(x); //*/

  if (k0<3.*k) { 
   _4D+=
   make([&](double pq, double pqd) { return make([&](double r, double rd) {
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double qm = km-q, qp = kp-q;
      double pm = km-p, pp = kp-p;
      double h = fabs(2.*km);
      //double h = fabs(2.*(k-pq));

      //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
      //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
      //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
      //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

      temp += F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qp/(q*q) );
      temp += F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pp/(p*p) );
      temp += F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qp/(q*q) );
      temp += F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pp/(p*p) );
      //double L = lga (qm*qp*pm*pp/SQR(p*q));
      //temp += F_14(p,k0-p)*(.5+f(q,s2))*L;
      //temp += F_14(q,k0-q)*(.5+f(p,s2))*L;
      //temp += F_25(p,k0-p)*(.5+f(q,s1))*L;
      //temp += F_25(q,k0-q)*(.5+f(p,s1))*L;
      return .5*temp/r;
  },  -km,km  )(y); }, 0.,k-km    )(x); 
   _4D+=
   make([&](double pq, double pqd) { return make([&](double r, double rd) {
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double qm = km-q, qp = kp-q;
      double pm = km-p, pp = kp-p;
      double h = fabs(k-pq-km);
      //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
      //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
      //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
      //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

      temp += lga(qm*pm/(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 2.*lga( q*(k0-q)/(qp*qm) )*F_14(p,k0-p)*(.5+f(r,s3)); // virtual
      temp += 2.*lga( p*(k0-p)/(pp*pm) )*F_14(q,k0-q)*(.5+f(r,s3)); // virtual
      return .5*temp/r;
  },  km,k-pq  )(y); }, 0.,k-km    )(x);
   _4D+=
   make([&](double pq, double pqd) { return make([&](double r, double rd) {
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double qm = km-q, qp = kp-q;
      double pm = km-p, pp = kp-p;
      double h = fabs(k-pq-km);
      //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
      //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
      //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
      //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

      temp += lga(qm*qp*pm*pp/SQR(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
      temp += 2.*lga( q*(k0-q)/(qp*qm) )*F_14(p,k0-p)*(.5+f(r,s3)); // virtual
      temp += 2.*lga( p*(k0-p)/(pp*pm) )*F_14(q,k0-q)*(.5+f(r,s3)); // virtual
      return .5*temp/r;
  },  pq-k,-km  )(y); }, 0.,k-km    )(x); }//*/

   res += _1A + _1B + _1D + _2A + _2B + _2C + _2D + _4A + _4B + _4D;
   //res += _4D;

  /*  #[ from '4' */

  } else {
    res += 0.;
  }
    // still don't have a good way to cater for NaNs
    if ( isinf(res)||isnan(res) ) { return 0.;}
    else { return fabs(K2)*res/k; }
    //else { return res/k; }
}
