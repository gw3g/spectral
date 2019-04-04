#include "core.hh"
#include "quad.hh"
#include "trapezoid.hh"
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
  };//*/
  double eval() {
    double res, err;
    //outer f1;
    //f1.f2.R = this;
    //integrate<outer> I(f1); // do the x-integral
    //res = go(I);
  double epsabs = 1e-3, epsrel = 1e-3;
size_t limit = 1e5;

  quad wsp1(limit);
  quad wsp2(limit);

  auto outer = make_gsl_function( [&](double x) {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
        return (this->integrand)(x,y);
             //return (this->integrand)(.5*(1.-x),y)+(this->integrand)(1.-.5*x,y);
        } );
    gsl_integration_qag(inner, .0+1e-10,1., epsabs, epsrel, 
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag(  outer, .0+1e-10,1., epsabs, epsrel, 
                         limit, 6, wsp2, &res, &err  );//*/

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
  res*= pow(k0,-_m)*pow(p,_m)*sgn(km);

  return res;
}

double rho11111::F_25(double q, double v) {
  double res;

  int s2 = (this->s)[2],
      s5 = (this->s)[5];

  int _n = this->n;

  double fq=f(q,s2), fv=f(v,s5);

  res = ( ((double) s2*s5)*exp(k0)-1. )*fq*fv;
  res*= pow(k0,-_n)*pow(q,_n)*sgn(km);

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
  double pp, pm, qp, qm;

  if (k0>k) {
  double _1A, _1B,      _1D,
         _2A, _2B, _2C, _2D, // 3 by reflection
         _4A, _4B,      _4D;

  //  #[ from '1'

  _1A =
  make([&](double q, double qd) {
       qm = km-q; qp = kp-q;
       return make([&](double p, double pd) {
          double temp = 0., r=k0-p-q;
          pm = km-p; pp = kp-p;
          if (k0>3.*k) {
            temp +=  lga(pp*qp/(pm*qm))*(
                                          F_123(p,q,r)+F_123(q,p,r)
                                        + F_345(p,q,r)+F_345(q,p,r)
                                        );
          };
          return temp/r;
  },  max(qp,q), km  )(y); },  k, km    )(x);//*/

  _1B =
  make([&](double q, double qd) {
       qm = km-q;
       return make([&](double p, double pd) {
          double temp = 0., r = k0-p-q;
          pm = km-p; 
          // p=km (!)
          //pm = (fabs(pm)<1e-1*( min(km,kp-q)-max(km-q,q) )) ? pd : pm;
          temp += lga( p*q/(pm*qm) )*(
                                       F_123(p,q,r)+F_123(q,p,r) 
                                     + F_345(p,q,r)+F_345(q,p,r) 
                                     );
          return temp/r;
  },  max(km-q,q), min(km,kp-q) )(y); }, 0.,min(km,.5*kp)   )(x);//*/

  _1D = 0.; // see 4B

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
  make([&](double p, double pd) {
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return make([&](double q, double qd) {
        double temp = 0., r = k0-p-q;
        qm = km-q; qp = kp-q;
        temp += lga(pm*qm/(pp*qp))*(
                                     F_123(p,q,r) + F_123(q,p,r)
                                   + F_345(p,q,r) + F_345(q,p,r)
                                   );
        //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
        return temp/r;
  },  km-p  )(y); },  kp )(x); //*/

  _2B =
  make([&](double p, double pd) { // [kp,+inf)
      pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return make([&](double q, double qd) { // [km-p,kp-p]
        double temp = 0., r = k0-p-q;
        qp = kp-q;
        temp += lga(p*q/(pp*qp))*(
                                   F_123(p,q,r) + F_123(q,p,r)
                                 + F_345(p,q,r) + F_345(q,p,r)
                                 );
        //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
        return temp/r;
  },  km-p,kp-p  )(y); },  kp )(x); //*/
  return _2A + _2B;

  _2C =
  make([&](double p, double pd) { // [km,kp]
      pm = km-p; pp = kp-p;
      //pm = (fabs(pm)<1e-1*k) ? pd : pm; // p=km (!)
      //pp = (fabs(pp)<1e-1*k) ? pd : pp; // p=kp (!)
      double l=k0-p;
      return make([&](double q, double qd) { // [km-p,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          double v=k0-q;
          //double q_ = (fabs(qm-p)<k0) ? qd : q; // q=0 (!)
          //temp += lga(pm*qm/(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
          //temp += lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          //
          temp += log( pm*qm/(p*q) )*(
                                       (.5+f(q,s2))*pow(q/k0,n)*F_14(p,k0-p)
                                     + (.5+f(q,s1))*pow(q/k0,m)*F_25(p,k0-p)
                                     );
          temp += log( pp*pm*q*v/( qm*qp*p*l ) )*(.5+f(r,s3))*(
                                       pow(q/k0,n)*F_14(p,k0-p)
                                     + pow(q/k0,m)*F_25(p,k0-p)
                                     );
          temp += log( l*v/(pp*qp) )*(
                                       (.5+f(v,s5))*pow(q/k0,n)*F_14(p,k0-p)
                                     + (.5+f(v,s4))*pow(q/k0,m)*F_25(p,k0-p)
                                     );
          return temp/r;
  },  km-p  )(y); },  km,kp )(x); //*/

  _2D =
  make([&](double p, double pd) { // [km,kp]
      return make([&](double q, double qd) { // [km-p,0]
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          double q_ = (fabs(q)<1e-1*fabs(km-p)) ? qd : q; // q=0 (!)
          //temp += lga( q_*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp += 2.*lga( q_*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( 
                                                             pow(q/k0,m)*F_14(p,k0-p)
                                                           + pow(q/k0,n)*F_25(p,k0-p)
                                                           );
          return temp/r;
  },  km-p, kp-p  )(y); },  km,kp    )(x); //*/

  //  #] from '2','3'\
      #[ from '4'

  _4A =
  make([&](double p, double pd) { // [kp,+inf)
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return make([&](double q, double qd) { // [p,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          //qp = (fabs(qp)<k0) ? qd : qp; // q=kp (!)
          temp += lga(pm*qm/(pp*qp))*(
                                       F_123(p,q,r) + F_123(q,p,r)
                                     + F_345(p,q,r) + F_345(q,p,r)
                                     );
          //return W_vi(p,q)*F_123(p,q,r)/r;
          return temp/r;
  },  p  )(y); },  kp    )(x); //*/

  _4B =
  make([&](double q, double qd) { // [km,kp]
      qm = km-q; qp = kp-q;
      //qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=km (!)
      //qp = (fabs(qp)<1e-1*k) ? qd : qp; //  =kp (!)
      double v = k0-q;
      return make([&](double p, double pd) { // [k0,+inf)
          double temp = 0., r = k0-p-q;
          pm = km-p; pp = kp-p;
          double l = k0-p;
          //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp += lga(pm*qm/(p*q))*(
                                   F_123(p,q,r) + F_123(q,p,r)
                                 + F_345(p,q,r) + F_345(q,p,r)
                                 );
          temp += 2.*lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*(
                                                         pow(p/k0,n)*F_14(q,k0-q) + 
                                                       + pow(p/k0,m)*F_25(q,k0-q) 
                                                       );
          //temp += log( pm*qm/(p*q) )*(
                                       //(.5+f(p,s2))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(p,s1))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( qp*qm*p*l/( pm*pp*q*v ) )*(.5+f(r,s3))*(
                                       //pow(p/k0,n)*F_14(q,k0-q)
                                     //+ pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( l*v/(pp*qp) )*(
                                       //(.5+f(l,s5))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(l,s4))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          return temp/r;
  },  k0  )(y); },  km,kp    )(x);

  if (k0>3.*k) {
  //
    _4B+=
    make([&](double q, double qd) { // [km,kp]
        qm = km-q; qp = kp-q;
        //qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=km (!)
        double v = k0-q;
        return make([&](double p, double pd) { // [k0+km-q,k0]
            double temp = 0., r = k0-p-q;
            pm = km-p; pp = kp-p;
            double l = k0-p; //l = (fabs(l)<1e-1*fabs(q-km)) ? pd : l; // p=k0 (!)
            temp += lga(pm*qm/(p*q))*( 
                F_123(p,q,r) + F_123(q,p,r) 
                );
            //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
            temp += lga( p*l/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
            temp += lga( p*l/(pp*pm) )*(.5+f(r,s3))*( F_14(k0-q,q) + F_25(q,k0-q) );
          //temp += log( pm*qm/(p*q) )*(
                                       //(.5+f(p,s2))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(p,s1))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( qp*qm*p*l/( pm*pp*q*v ) )*(.5+f(r,s3))*(
                                       //pow(p/k0,n)*F_14(q,k0-q)
                                     //+ pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( l*v/(pp*qp) )*(
                                       //(.5+f(l,s5))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(l,s4))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
            return temp/r;
    },  k0+km-q,k0  )(y); },  km,kp    )(x);
  //
    _4B+=
    make([&](double q, double qd) { // [km,kp]
        qm = km-q; qp = kp - q;
        //qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=km (!)
        double v = k0-q;
        return make([&](double p, double pd) { // [2kp-q,k0+km-q]
            double temp = 0., r = k0-p-q;
            pm = km-p; pp = kp-p;
            double l = k0-p; //l = (fabs(l)<1e-1*(km-k)) ? pd : l;
            temp += lga(pm*qm/(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
            temp -= lga( pm*qm/((k0-q)*l) )*( F_123(l,k0-q,-r) + F_123(k0-q,l,-r) ); // from 1D & 1C
            temp += lga( p*l/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
            temp += lga( p*l/(pp*pm) )*(.5+f(r,s3))*( F_14(k0-q,q) + F_25(q,k0-q) );
          //temp += log( pm*qm/(p*q) )*(
                                       //(.5+f(p,s2))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(p,s1))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( qp*qm*p*l/( pm*pp*q*v ) )*(.5+f(r,s3))*(
                                       //pow(p/k0,n)*F_14(q,k0-q)
                                     //+ pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( l*v/(pp*qp) )*(
                                       //(.5+f(l,s5))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(l,s4))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
            return temp/r;
    },  2.*kp-q,k0+km-q  )(y); },  km,kp    )(x);
  //
  } else if (k0<3.*k) {
  //
    _4B+=
    make([&](double p, double pd) { // [kp,k0]
        pm = km-p; pp = kp-p;
        //pp = (fabs(pp)<1e-1*fabs(k0-kp)) ? pd : pp; // p=kp (!)
        double l = k0-p; //l = (fabs(l)<1e-1*fabs(k0-kp)) ? pd : l; // p=k0 (!)
        return make([&](double q, double qd) { // [k0+km-p,kp]
            double temp = 0., r = k0-p-q;
            qm = km-q; qp = kp-p;
            double v = k0-q;
            temp += lga(pm*qm/(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
            temp += lga( p*l/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
            temp += lga( p*l/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          //temp += log( pm*qm/(p*q) )*(
                                       //(.5+f(p,s2))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(p,s1))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( qp*qm*p*l/( pm*pp*q*v ) )*(.5+f(r,s3))*(
                                       //pow(p/k0,n)*F_14(q,k0-q)
                                     //+ pow(p/k0,m)*F_25(q,k0-q)
                                     //);
          //temp += log( l*v/(pp*qp) )*(
                                       //(.5+f(l,s5))*pow(p/k0,n)*F_14(q,k0-q)
                                     //+ (.5+f(l,s4))*pow(p/k0,m)*F_25(q,k0-q)
                                     //);
            return temp/r;
    },  k0+km-p, kp  )(y); },  kp,k0    )(x); //*/
  }
  _4B*=2.;

  _4D =
  make([&](double pq, double pqd) { // [max(0,k-km),k]
      return make([&](double r, double rd) { // [p-q-k,k-p+q]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          pm = km-p; qm = km-q;
          pp = kp-p; qp = kp-q;
          double l=k0-p, p2=p*p, l2=l*l;
          double v=k0-q, q2=q*q, v2=v*v;

          double h = fabs(k-pq);
          qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;

          temp += pow(q,n)*F_14(p,l)*( (.5+f(q,s2))*lga( qm*qp/q2 ) - (.5+f(v,s5))*lga( qm*qp/v2 ) ); // orig
          temp += pow(p,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*pp/p2 ) - (.5+f(l,s4))*lga( pm*pp/l2 ) );

          temp += pow(p,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*pp/p2 ) - (.5+f(l,s5))*lga( pm*pp/l2 ) ); // p <--> q, l <--> v
          temp += pow(q,m)*F_25(p,l)*( (.5+f(q,s1))*lga( qm*qp/q2 ) - (.5+f(v,s4))*lga( qm*qp/v2 ) );

          temp -= pow(v,n)*F_14(l,p)*( (.5+f(v,s2))*lga( qm*qp/v2 ) - (.5+f(q,s5))*lga( qm*qp/q2 ) ); // p <--> l, q <--> v
          temp -= pow(l,m)*F_25(v,q)*( (.5+f(l,s1))*lga( pm*pp/l2 ) - (.5+f(p,s4))*lga( pm*pp/p2 ) );

          temp -= pow(l,n)*F_14(v,q)*( (.5+f(l,s2))*lga( pm*pp/l2 ) - (.5+f(p,s5))*lga( pm*pp/p2 ) ); // p <--> v, q <--> l
          temp -= pow(v,m)*F_25(l,p)*( (.5+f(v,s1))*lga( qm*qp/v2 ) - (.5+f(q,s4))*lga( qm*qp/q2 ) );


          //temp += F_14(q,v)*(.5+f(p,s2))*lga( pm*pp/p2 )*pow(p,n);
          //temp += F_25(p,l)*(.5+f(q,s1))*lga( qm*qp/q2 )*pow(q,m);
          //temp += F_25(q,v)*(.5+f(p,s1))*lga( pm*pp/p2 )*pow(p,m);

          //temp -= F_14(v,q)*(.5+f(l,s2))*lga( pm*pp/l2 )*pow(l,n);
          //temp -= F_14(l,p)*(.5+f(v,s2))*lga( qm*qp/v2 )*pow(v,n);
          //temp -= F_25(v,q)*(.5+f(l,s1))*lga( pm*pp/l2 )*pow(l,m);
          //temp -= F_25(l,p)*(.5+f(v,s1))*lga( qm*qp/v2 )*pow(v,m);
          return 2.*.5*temp/r;
  },  0.,k-pq  )(y); }, 0.,k    )(x); //*/

  _1D += // corner, 1D+C
  make([&](double pq, double pqd) { // [k-min(k,km),k+min(k,km)]
      return make([&](double r, double rd) { // [|pq-k|,min(k,km)]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          pm = km-p; qm = km-q;
          pp = kp-p; qp = kp-q;
          double h = fabs( min(k,km)-fabs(pq-k) );
          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
          //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

          temp += lga( pp*qp/(q*p) )*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( pp*qp/((k0-p)*(k0-q)))*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );
          temp += 2.*lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return .5*temp/r;
  },  fabs(pq-k),min(k,km)  )(y); }, k-min(k,km),k+min(k,km)  )(x); //*/

  if (k0<3.*k) { 
  //
    _4D+=
    make([&](double pq, double pqd) { // [0,k-km]
        return make([&](double r, double rd) { // [-km,km]
            double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
            pm = km-p; qm = km-q;
            pp = kp-p; qp = kp-q;
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
            return 2.*.5*temp/r;
    },  -km,km  )(y); }, 0.,k-km    )(x); 
  //
    _4D+=
    make([&](double pq, double pqd) { // [0,k-km]
        return make([&](double r, double rd) { // [km,k-pq]
            double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
            pm = km-p; qm = km-q;
            pp = kp-p; qp = kp-q;
            double h = fabs(k-pq-km);
            //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
            //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
            //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
            //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

            temp += lga(qm*pm/(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
            temp += lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
            temp += lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
            return 2.*.5*temp/r;
    },  km,k-pq  )(y); }, 0.,k-km    )(x);
  //
    _4D+=
    make([&](double pq, double pqd) { // [0,k-km]
        return make([&](double r, double rd) { // [pq-k,-km]
            double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
            pm = km-p; qm = km-q;
            pp = kp-p; qp = kp-q;
            double h = fabs(k-pq-km);
            //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
            //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
            //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
            //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;

            temp += lga(qm*qp*pm*pp/SQR(p*q))*( F_123(p,q,r) + F_123(q,p,r) );
            temp += lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
            temp += lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
            return .5*temp/r;
    },  pq-k,-km  )(y); }, 0.,k-km    )(x); //*/
  }

  res += _1A + _1B + _2A + _2B + _2C + _4A ;
  //res += _1A + _1B + _1D + _2A + _2B + _2C + _2D + _4A + _4B + _4D;

  } else {

  double _2A,
         _4A, _4B,      _4D,
         _5A, _5B,
         _6A, _6B, _6C; // 7 by reflection

  //  #[ from '2'

  _2A =
  make([&](double p, double pd) { // [km,0]
      pm = km-p; pp = kp-p;
      return make([&](double q, double qd) { // [km,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          //temp += lga( pp*pm/(p*(k0-p)) )*( F_123(p,q,r) + F_123(q,p,r) );
          // ^ ZERO?
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  km  )(y); },  km,0.    )(x);

  _2A +=
  make([&](double p, double pd) { // [0,k0]
      pm = km-p, pp = kp-p;
      return make([&](double q, double qd) { // [km,+inf)
          double temp = 0., r = k0-p-q;
          //qm = km-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          //temp += lga( p*q/(pm*qm) )*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  km  )(y); },  0.,k0    )(x);

  _2A +=
  make([&](double p, double pd) { // [k0,kp]
      pm = km-p, pp = kp-p;
      return make([&](double q, double qd) { // [k0+km-p,+inf)
      double temp = 0., r = k0-p-q;
      temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
      return temp/r;
  },  k0+km-p  )(y); },  k0,kp    )(x); //*/

  //  #] from '2'\
      #[ from '4'

  _4A =
  make([&](double p, double pd) { // [kp,+inf)
      return make([&](double q, double qd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          temp += W_vi(p,q)*( F_123(p,q,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp  )(y); },  kp    )(x); //*/

  _4B =
  make([&](double q, double qd) { // [k0,kp]
      qm = km-q; qp = kp-q;
      return make([&](double p, double pd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp  )(y); },  k0,kp    )(x); //*/

  _4B +=
  make([&](double q, double qd) { // [0,k0]
      qm = km-q; qp = kp-q;
      return make([&](double p, double pd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp  )(y); },  0.,k0    )(x);

  _4B +=
  make([&](double q, double qd) { // [km,0]
      qm = km-q; qp = kp-q;
      return make([&](double p, double pd) { // [kp-q,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp-q  )(y); }, km, 0.  )(x); //*/

  _4D =
  make([&](double q, double qd) { // [0,kp]
      qp = kp-q;
      qp = (fabs(qp)<1e-1*kp) ? qd : qp; // q=kp (!)
      return make([&](double p, double pd) { // [kp-q,kp]
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          pp = (fabs(pp)<1e-1*q) ? pd : pp; // p=kp (!)
          temp += lga( pp*qp/(p*q) )*( F_123(p,q,r) );
          return temp/r;
  },  kp-q,kp  )(y); },  0.,kp    )(x); //*/

  //  #] from '4'\
      #[ from '5'

  _5A =
  make([&](double q, double qd) { // (-inf,km]
      qm = km-q; qp = kp-q;
      return make([&](double p, double pd) { // (-inf,km]
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  km  )(y); },  km    )(x); //*/

  _5B =
  make([&](double q, double qd) { // [km,0]
      qm = km-q;
      qm = (fabs(qm)<1e-1*fabs(km)) ? qd : qm; // q=km (!)
      return make([&](double p, double pd) { // [km,km-q]
          double temp = 0., r = k0-p-q;
          pm = km-p;
          pm = (fabs(pm)<1e-1*fabs(q)) ? pd : pm; // p=kp (!)
          temp += lga( pm*qm/(p*q) )*( F_123(p,q,r) );
          return temp/r;
  },  km,km-q  )(y); },  km,0    )(x); //*/

  //  #] from '5'\
      #[ from '6'

  _6A =
  make([&](double p, double pd) { // [k,+inf)
      pm = km-p, pp = kp-p;
      return make([&](double q, double qd) { // [kp-p,km]
          double temp = 0., r = k0-p-q;
          qm = km-q, qp = kp-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  kp-p,km  )(y); },  k    )(x); //*/

  _6B =
  make([&](double pq, double pqd) { // [k-km,+inf)
      return make([&](double r, double rd) { // [km,-km]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          qm = km-q, qp = kp-q;
          pm = km-p, pp = kp-p;
          double h = 2.*fabs(km);

          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
          //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
          //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          //temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          //temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );

          temp -= F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qm/(q*q) );
          temp -= F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pm/(p*p) );
          temp -= F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qm/(q*q) );
          temp -= F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pm/(p*p) );

          return .5*temp/r;
  },  km,-km  )(y); }, k-km    )(x); //*/

  _6B +=
  make([&](double p, double pd) { // [kp,+inf)
      pm = km-p, pp = kp-p;
      return make([&](double q, double qd) { // [km-p,km-p+k0]
          double temp = 0., r = k0-p-q;
          qm = km-q, qp = kp-q;
          double h = fabs( k0 );

          //pp = (fabs(pp)<1e-1*h) ? pd : pp;
          //qm = (fabs(qm)<1e-1*h) ? qd : qm;
          //pm = (fabs(pm)<1e-1*h) ? pd : pm;
          //qp = (fabs(qp)<1e-1*h) ? qd : qp;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );

          //temp += F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qp/(q*q) );
          //temp += F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pp/(p*p) );
          //temp += F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qp/(q*q) );
          //temp += F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pp/(p*p) );
          return temp/r;
  },  km-p,k0+km-p  )(y); },  kp    )(x); //*/

  _6B +=
  make([&](double pq, double pqd) { // [k,k-km]
      return make([&](double r, double rd) { // [k-pq,pq-k]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          qm = km-q, qp = kp-q;
          pm = km-p, pp = kp-p;
          double h = 2.*fabs( pq-k );

          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
          //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
          //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          //temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          //temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );

          temp -= F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qm/(q*q) );
          temp -= F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pm/(p*p) );
          temp -= F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qm/(q*q) );
          temp -= F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pm/(p*p) );
          return .5*temp/r;
  },  k-pq,pq-k  )(y); }, k,k-km    )(x); //*/

  _6B +=
  make([&](double pq, double pqd) { // [k+km,k-km]
      return make([&](double r, double rd) { // [|p-q-k|,-km]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          qm = km-q, qp = kp-q;
          pm = km-p, pp = kp-p;
          double h = fabs( km+fabs(pq-k) );

          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
          //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
          //temp -= W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          //temp += W_vi(k0-q,k0-p)*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );
          temp += lga( pm*qm/(q*p) )*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( pm*qm/((k0-p)*(k0-q)))*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );
          //temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );

          //temp += lga( pp*qp/(q*p) )*( F_123(p,q,r) + F_123(q,p,r) );
          //temp -= lga( pp*qp/((k0-p)*(k0-q)))*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );
          //temp -= 2.*lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp += lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );

          return .5*temp/r;
  },  fabs(pq-k),-km  )(y); }, k+km,k-km  )(x); //*/

  _6C =
  make([&](double q, double qd) { // (-inf,-k]
      qm = km-q, qp = kp-q;
      return make([&](double p, double pd) { // [kp,km-q]
          double temp = 0., r = k0-p-q;
          pm = km-p, pp = kp-p;
          double h = fabs( km-q-kp );

          //pp = (fabs(pp)<1e-1*h) ? pd : pp;
          //qm = (fabs(qm)<1e-1*h) ? qd : qm;
          //pm = (fabs(pm)<1e-1*h) ? pd : pm;
          //qp = (fabs(qp)<1e-1*h) ? qd : qp;

          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  kp,km-q  )(y); },  -k    )(x); //*/

  //  #] from '6'

    //res += _6B;
    res += _2A + _4A + _4B + _4D + _5A + _5B + _6A + _6B + _6C;
  }
    // still don't have a good way to cater for NaNs
    if ( isinf(res)||isnan(res) ) { return 0.;}
    else { return fabs(K2)*res/k; }
    //else { return res/k; }
}
