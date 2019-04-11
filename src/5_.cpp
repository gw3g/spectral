#include "core.hh"
#include "quad.hh"
#include "map.hh"

/*--------------------------------------------------------------------*/

struct rho11110 : Master {
  double integrand(double,double); // supported on [0,1]x[0,1]

  double F_123(double,double,double); // Cut: f1f2f3/f0   }- real
  double F_14(double,double);         //        f1f4/f0   }- virtual


  /*struct inner {
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
  };//*/
  double eval() {
    double res, err;
    //outer f1;
    //f1.f2.R = this;
   //integrate<outer> I(f1); // do the x-integral
    //res = go(I) + ( (k0>k) ? -K2/8. : 0. );
  double epsabs = 1e-2, epsrel = 0;
  size_t limit = 1e5;

  quad wsp1(limit);
  quad wsp2(limit);

  auto outer = make_gsl_function( [&](double x) {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
        return (this->integrand)(x,y);
        } );
    gsl_integration_qag(inner, .0+1e-10,1., epsabs, epsrel, 
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag(  outer, .0+1e-10,1., epsabs, epsrel, 
                         limit, 6, wsp2, &res, &err  );

    //return (( res + ( (k0>k) ? -K2/8. : 0. ) ))*CUBE(OOFP);//*/
    return (( res ))*CUBE(OOFP);
  }
  rho11110(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=5; }
};
// function for MAIN
Master* _11110(int m, int n, int s[3]) {
  Master *R =  new rho11110(m,n,s); return R;
}

/*#include <fstream>
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
}//*/

/*--------------------------------------------------------------------*/
// all the thermal weights for cuts from this topology

double rho11110::F_123(double p, double q, double r) { 
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

double rho11110::F_14(double p, double l) {
  double res;

  int s1 = (this->s)[1],
      s4 = (this->s)[4];

  int _m = this->m;

  double fp=f(p,s1), fl=f(l,s4);

  res = ( ((double) s1*s4)*exp(k0)-1. )*fp*fl;
  res*= pow(k0,-_m)*pow(p,_m)*sgn(km);

  return res;
}



/*--------------------------------------------------------------------*/

double rho11110::integrand(double x, double y) 
{
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
  remap([&](double p, double pd) { // [.5kp,km]
       pm = km-p; pp = kp-p;
       pm = (fabs(pm)<1e-1*fabs(km-.5*kp)) ? pd : pm; // p=km (!)
       return remap([&](double q, double qd) { // [kp-p,p]
          double temp = 0., r=k0-p-q;
          qm = km-q; qp = kp-q;
          double l=k0-p, v=k0-q;
          if (k0>3.*k) {
            temp +=  lga( kp*pp/(km*pm) )*F_123(p,q,r)/l ;
            temp +=  lga( kp*qp/(km*qm) )*F_123(q,p,r)/v ; // 1A
          } else {
            // this si to give correct log (km r/(kp-p)q),
            // when combined w/ region 4D.
            temp +=  lga( pp*r/(km*q) )*F_123(p,q,r)/l ;
            temp +=  lga( qp*r/(km*p) )*F_123(q,p,r)/v ; // 1A'
            temp *= -1.; // (signed integration region)
            // p -> -p and q -> k0-q (missing part of 3D)
            temp +=  lga( km*kp/((kp+p)*(km+p)) )*F_123(-p,v,p+q)/(k0+p) ;
            temp +=  lga( km*kp/((kp+q)*(km+q)) )*F_123(-q,l,p+q)/(k0+q) ;
          }
          return temp;
  },  kp-p,p  )(y); },  .5*kp, km    )(x);//*/
  //\
  return _1A;
  _1B =
  remap([&](double p, double pd) { // [.5km,km]
       pm = km-p; 
       pm = (fabs(pm)<1e-1*( .5*km )) ? pd : pm; // p=km (!)
       return remap([&](double q, double qd) { // [km-p,min(p,kp-p)]
          double temp = 0., r = k0-p-q;
          qm = km-q;
          double l=k0-p, v=k0-q;
          temp += lga( kp*q/(r*pm) )*F_123(p,q,r)/l ;
          temp += lga( kp*p/(r*qm) )*F_123(q,p,r)/v ;
          return temp;
  },  km-p,min(p,kp-p) )(y); }, .5*km,km   )(x);//*/
  //\
  return _1B;
  //
  _1D = // corner, 1D+C
  remap([&](double pq, double pqd) { // [k-min(k,km),k+min(k,km)]
      return remap([&](double r, double rd) { // [|pq-k|,min(k,km)]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          pm = km-p; qm = km-q;
          pp = kp-p; qp = kp-q;
          double l=k0-p, v=k0-q;
          //\
          Padding for log singularities:
          double h = fabs( min(k,km)-fabs(pq-k) );
          if (pq>k) { pp = (fabs(pp)<1e-1*h) ? rd/2. : pp; }
          if (pq<k) { qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; }
          
          temp += lga( pp*r/(km*q) )*F_123(p,q,r)/l ; // 1D
          temp += lga( kp*r/(p*qm) )*F_123(q,p,r)/v ; // 1C
          return .5*temp;
  },  fabs(pq-k),min(k,km)  )(y); }, k-min(k,km),k+min(k,km)  )(x); //*/\
  return _1D;

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
  remap([&](double p, double pd) { // [kp,k0+km]
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // (-inf,km-p]
        double temp = 0., r = k0-p-q;
        double l=k0-p, v=k0-q;
        qm = km-q; qp = kp-q;
        temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ;
        temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ;
        return temp;
  },  km-p  )(y); },  kp,k0+km )(x); //*/
  _2A +=
  remap([&](double p, double pd) { // [k0+km,+inf)
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // (-inf,km-p]
        double temp = 0., r = k0-p-q;
        double l=k0-p, v=k0-q;
        qm = km-q; qp = kp-q;
        temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ;
        temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ;
        return temp;
  },  km-p  )(y); },  k0+km )(x); //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) { // [kp,k0+km]
      pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // [km-p,kp-p]
        double temp = 0., r = k0-p-q;
        qp = kp-q;
        double l=k0-p, v=k0-q;
        temp += lga( km*q/(pp*r) )*F_123(p,q,r)/l ;
        temp += lga( km*p/(qp*r) )*F_123(q,p,r)/v ;
        return temp;
  },  km-p,kp-p  )(y); },  kp,k0+km )(x); //*/
  _2B +=
  remap([&](double p, double pd) { // [k0+km,+inf)
      pp = kp-p;
      return remap([&](double q, double qd) { // [km-p,kp-p]
        double temp = 0., r = k0-p-q;
        qp = kp-q;
        double l=k0-p, v=k0-q;
        temp += lga( km*q/(pp*r) )*F_123(p,q,r)/l ;
        temp += lga( km*p/(qp*r) )*F_123(q,p,r)/v ;
        return temp;
  },  km-p,kp-p  )(y); },  kp,k0+km )(x); //*/
  //\
  return _2B;
  _2C =
  remap([&](double p, double pd) { // [km,kp]
      pm = km-p; pp = kp-p;
      pm = (fabs(pm)<1e-1*k) ? pd : pm; // p=km (!)
      pp = (fabs(pp)<1e-1*k) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // (-inf,km-p]
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          double l=k0-p, v=k0-q;
          temp += lga( pm*r/(kp*q) )*F_123(p,q,r)/l ;
          temp += lga( km*r/(qp*p) )*F_123(q,p,r)/v ;
          return temp;
  },  km-p  )(y); },  km,kp )(x); //*/
  //\
  return _2C;

  //  #] from '2','3'\
      #[ from '4'

  _4A =
  remap([&](double q, double qd) { // [kp,k0]
      return remap([&](double p, double pd) { // [2k0-q,+inf)
          double temp = 0., r = k0-p-q;
          pm = km-p; pp = kp-p;
          qm = km-q; qp = kp-q;
          double l=k0-p, v=k0-q;
          //qp = (fabs(qp)<k0) ? qd : qp; // q=kp (!)
          temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ;
          temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ;
          return temp;
  },  2.*k0-q  )(y); }, kp,k0    )(x); //*/
  _4A +=
  remap([&](double q, double qd) { // [kp,k0]
      return remap([&](double p, double pd) { // [q,2kp-q]
          double temp = 0., r = k0-p-q;
          pm = km-p; pp = kp-p;
          qm = km-q; qp = kp-q;
          double l=k0-p, v=k0-q;
          //qp = (fabs(qp)<k0) ? qd : qp; // q=kp (!)
          temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ;
          temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ;
          return temp;
  },  q,2.*k0-q  )(y); },  kp,k0    )(x); //*/
  _4A +=
  remap([&](double q, double qd) { // [k0,+inf)
      return remap([&](double p, double pd) { // [q,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          pm = km-p; pp = kp-p;
      //pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
          double l=k0-p, v=k0-q;
          //qp = (fabs(qp)<k0) ? qd : qp; // q=kp (!)
          temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ;
          temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ;
          return temp;
  },  q  )(y); },  k0    )(x); //*/
  //
  return _4A;

    res += _1A + _1B;
    //res += _1A + _1B + _1D + _2A + _2B + _2C + _2D + _4A + _4B + _4D;

  } else {

    //res += _2A + _4A + _4B + _4D + _5A + _5B + _6A + _6B + _6C;

  }

    // still don't have a good way to cater for NaNs
  if ( isinf(res)||isnan(res) ) { return 0.;}
  else { return res/(k); }
}
