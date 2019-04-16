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
    // check included cases:
    if      ( (m==0)&&(n==0) ) {  }
    else if ( (m==1)&&(n==0) ) {  }
    else if ( (m==0)&&(n==1) ) {  }
    else                       { return 0.;}

    double res, err;
    //outer f1;
    //f1.f2.R = this;
   //integrate<outer> I(f1); // do the x-integral
    //res = go(I) + ( (k0>k) ? -K2/8. : 0. );
    double epsabs = 1e-2, epsrel = 1e-2;
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
  rho11110 *R =  new rho11110(m,n,s);
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

  double _Virt,
         _1A, _1B,      _1D,
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
            // this is to give correct log[km r/(kp-p)q],
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
          temp += lga( kp*p/(r*qm) )*F_123(q,p,r)/v ; // 1B
          return temp;
  },  km-p,min(p,kp-p) )(y); }, .5*km,km   )(x);//*/
  //\
  return _1B;
  //
  _1D = 
  remap([&](double p, double pd) { // [max(k,km),kp]
      return remap([&](double q, double qd) { // [kp-p,km]
          double temp = 0., r=k0-p-q;
          qm = km-q;
          pp = kp-p;
          double l=k0-p, v=k0-q;
          //\
          Padding for log singularities:
          double h=fabs(km-kp+p);
          qm = (fabs(qm)<1e-1*h) ? qd : qm;

          temp += lga( pp*r/(km*q) )*F_123(p,q,r)/l ; // 1D
          temp += lga( kp*r/(p*qm) )*F_123(q,p,r)/v ; // 1C
          return temp;
  },  kp-p,km  )(y); }, max(k,km),kp  )(x); //*/
  //\
  return _1D;

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
  remap([&](double p, double pd) { // [kp,k0+km]
      // the sing. @ p=k0 is at the midpoint of the interval
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<1e-1*(k0-k)) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // (-inf,km-p]
        double temp = 0., r = k0-p-q;
        double l=k0-p, v=k0-q;
        qm = km-q; qp = kp-q;
        temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ; // 2A
        temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ; // 3A
        return temp;
  },  km-p  )(y); },  kp,k0+km )(x); //*/
  _2A +=
  remap([&](double p, double pd) { // [k0+km,+inf)
      pm = km-p; pp = kp-p;
      return remap([&](double q, double qd) { // (-inf,km-p]
        double temp = 0., r = k0-p-q;
        double l=k0-p, v=k0-q;
        qm = km-q; qp = kp-q;
        temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ; // 2A
        temp += lga( km*qm/(kp*qp) )*F_123(q,p,r)/v ; // 3A
        return temp;
  },  km-p  )(y); },  k0+km )(x); //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) { // [kp,k0+km]
      // the sing. @ p=k0 is at the midpoint of the interval
      pp = kp-p;
      pp = (fabs(pp)<1e-1*(k0-k)) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // [km-p,kp-p]
        double temp = 0., r = k0-p-q;
        qp = kp-q;
        double l=k0-p, v=k0-q;
        temp += lga( km*q/(pp*r) )*F_123(p,q,r)/l ; // 2B
        temp += lga( km*p/(qp*r) )*F_123(q,p,r)/v ; // 3B
        return temp;
  },  km-p,kp-p  )(y); },  kp,k0+km )(x); //*/
  _2B +=
  remap([&](double p, double pd) { // [k0+km,+inf)
      pp = kp-p;
      return remap([&](double q, double qd) { // [km-p,kp-p]
        double temp = 0., r = k0-p-q;
        qp = kp-q;
        double l=k0-p, v=k0-q;
        temp += lga( km*q/(pp*r) )*F_123(p,q,r)/l ; // 2B
        temp += lga( km*p/(qp*r) )*F_123(q,p,r)/v ; // 3B
        return temp;
  },  km-p,kp-p  )(y); },  k0+km )(x); //*/
  //\
  return _2B;
  _2C =
  remap([&](double p, double pd) { // [km,kp]
      pm = km-p;
      pm = (fabs(pm)<1e-1*k) ? pd : pm; // p=km (!)
      return remap([&](double q, double qd) { // (-inf,km-p]
          double temp = 0., r = k0-p-q;
          qp = kp-q;
          double l=k0-p, v=k0-q;
          temp += lga( pm*r/(kp*q) )*F_123(p,q,r)/l ; // 2B
          temp += lga( km*r/(qp*p) )*F_123(q,p,r)/v ; // 3B
          return temp;
  },  km-p  )(y); },  km,kp )(x); //*/
  //\
  return _2C;
  _2D =
  remap([&](double q, double qd) { // [-min(k,km),0]
      qm = km-q; qp = kp-q;
      return remap([&](double p, double pd) { // [km-q,kp]
          double temp = 0., r = k0-p-q;
          double l=k0-p, v=k0-q;
          temp += lga( km*kp/(qp*qm) )*F_123(q,p,r)/v ;// 3D (+part from 1A)
          temp += lga( km*kp/((kp+q)*(km+q)) )*F_123(-q,l,p+q)/(k0+q) ; // 1D'
          return temp;    // 2D+1C' are zero
  },  km-q,kp  )(y); },  -min(k,km),0. )(x); //*/
  //\
  return _2D;

  //  #] from '2','3'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) { // [kp,k0+km]
      // the sing. @ p=k0 is at the midpoint of the interval
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<1e-1*(k0-k)) ? pd : pp; // q=kp (!)
      return remap([&](double q, double qd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          double l=k0-p;
          temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ; // no reflection
          return temp;
  },  kp  )(y); }, kp,k0+km    )(x); //*/
  _4A +=
  remap([&](double p, double pd) { // [k0+km,+inf)
      pm = km-p; pp = kp-p;
      return remap([&](double q, double qd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          double l=k0-p;
          temp += lga( km*pm/(kp*pp) )*F_123(p,q,r)/l ; // no reflection
          return temp;
  },  kp  )(y); },  k0+km    )(x); //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) { // [km,kp]
      qm = km-q;
      qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=kp (!)
      return remap([&](double p, double pd) { // [kp,k0+km]
          // the sing. @ p=k0 is at the midpoint of the interval
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          double l=k0-p, v=k0-q;
          temp += lga( km*r/(q*pp) )*F_123(p,q,r)/l ; // 4B
          temp += lga( qm*r/(p*kp) )*F_123(q,p,r)/v ; // 4C
          return temp;
  },  kp,k0+km  )(y); }, km,kp    )(x); //*/
  _4B +=
  remap([&](double q, double qd) { // [km,kp]
      qm = km-q;
      qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=kp (!)
      return remap([&](double p, double pd) { // [kp,k0+km]
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          double l=k0-p, v=k0-q;
          temp += lga( km*r/(q*pp) )*F_123(p,q,r)/l ; // 4B
          temp += lga( qm*r/(p*kp) )*F_123(q,p,r)/v ; // 4C
          return temp;
  },  k0+km  )(y); }, km,kp    )(x); //*/
  //\
  return _4B;
  _4D =
  remap([&](double pq, double pqd) { // [0,k]
      return remap([&](double r, double rd) { // [-k+p-q,k-p+q]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          double l=k0-p, v=k0-q;

          temp += lga( r/q )*F_123(p,q,r)/l ; // 4D+1B+1B'
          temp += lga( r/p )*F_123(q,p,r)/v ; // (reflected)
          return temp; // no .5, from log(...^2)
  },  -k+pq,k-pq )(y); }, 0.,k    )(x); //*/
  //\
  return _4D;

  //  #] from '4',\
      #[ from virtual corrections

  _Virt =
  remap([&](double p, double pd) { // [km,kp]
      double tempp = 0.;
      double l=k0-p;
      pp = kp-p; pm = km-p;
      pp = (fabs(pp)<1e-1*k) ? pd : pp; // q=kp (!)
      pm = (fabs(pm)<1e-1*k) ? pd : pm; // q=km (!)
      tempp += remap([&](double q, double qd) { // [0,2l]
          double temp = 0., r=k0-p-q, rp=k0-p+q;

          temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*pow(q/k0,n) ;
          temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*pow(-q/k0,n) ;
          return temp;
      },  0.,2.*(k0-p) )(y); 
      tempp += remap([&](double q, double qd) { // [2l,+inf)
          double temp = 0., r=k0-p-q, rp=k0-p+q;

          temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*pow(q/k0,n) ;
          temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*pow(-q/k0,n) ;
          return temp;
      },  2.*(k0-p) )(y); 

      return F_14(p,l)*(
                      ( lga(4.*pm*pp/(k*k))+lga(l*l/K2)
                    + ( (m+n==1) ? 3.5 : 3. ) // 3.5 if m=1 !!!!
                    )*pow( .5*l/k0, n )
                    + 2.*s3*tempp/l
                    );
  }, km,kp    )(x); //*/
  //\
  return .25*_Virt/k;

    //res += _Virt;
    res += _1A + _1B + _1D + _2A + _2B + _2C + _2D + _4A + _4B + _4D + _Virt;

  } else {

  double _Virt,
         _2A, _2B,
         _4A, _4B,      _4D,
         _5A, _5B,
         _6A, _6B, _6C; // 7 by reflection

  //  #[ from '2'

  _2A =
  remap([&](double p, double pd) { // [km,0]
      pm = km-p; pp = kp-p;
      return remap([&](double q, double qd) { // (-inf,km]
          double temp = 0., r = k0-p-q;
          double l=k0-p;
          temp += lga( pp*pm/(kp*km) )*F_123(p,q,r)/l ;   // 2A
          return temp;                                    // 3A = 0
  },  km  )(y); },  km,0.    )(x);
  _2A +=
  remap([&](double p, double pd) { // [0,k0+km]
      pm = km-p, pp = kp-p;
      return remap([&](double q, double qd) { // (-inf,km-p]
          double temp = 0., r = k0-p-q;
          double l=k0-p;
          temp += lga( pp*pm/(kp*km) )*F_123(p,q,r)/l ;   // 2A
          return temp;                                    // 3A = 0
  }, km-p  )(y); },  0.,k0+km    )(x);
  _2A +=
  remap([&](double p, double pd) { // [k0+km,kp]
      pm = km-p, pp = kp-p;
      return remap([&](double q, double qd) { // (-inf,km-p]
          double temp = 0., r = k0-p-q;
          double l=k0-p;
          temp += lga( pp*pm/(kp*km) )*F_123(p,q,r)/l ;   // 2A
          return temp;                                    // 3A = 0
  },  km-p  )(y); },  k0+km,kp    )(x); //*/\
  return _2A;
  _2B  =
  remap([&](double p, double pd) { // [0,k0+km]
      pp = kp-p;
      return remap([&](double q, double qd) { // [km-p,km]
          double temp = 0., r = k0-p-q;
          qm = km-q;
          qm = (fabs(qm)<1e-1*p) ? qd : qm;
          double l=k0-p, v=k0-q;
          temp += lga( q*pp/(r*km) )*F_123(p,q,r)/l ;     // 2B
          temp += lga( p*kp/(r*qm) )*F_123(q,p,r)/v ;     // 3B
          return temp;
  },  km-p,km  )(y); },  0.,k0+km    )(x);
  _2B +=
  remap([&](double p, double pd) { // [k0+km,kp]
      pp = kp-p;
      return remap([&](double q, double qd) { // [km-p,km]
          double temp = 0., r = k0-p-q;
          qm = km-q;
          qm = (fabs(qm)<1e-1*p) ? qd : qm;
          double l=k0-p, v=k0-q;
          temp += lga( q*pp/(r*km) )*F_123(p,q,r)/l ;     // 2B
          temp += lga( p*kp/(r*qm) )*F_123(q,p,r)/v ;     // 3B
          return temp;
  },   km-p,km  )(y); },  k0+km,kp    )(x); //*/\
  return _2B;

  //  #] from '2'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) { // [kp,+inf)
      pm = km-p;
      return remap([&](double q, double qd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          double l=k0-p;
          temp += lga( q*pm/(r*kp) )*F_123(p,q,r)/l ;     // 4A
          return temp;
  },  kp  )(y); },  kp    )(x); //*/\
  return _4A;
  _4B =
  remap([&](double q, double qd) { // [km,kp]
      qp = kp-q; qm = km-q;
      return remap([&](double p, double pd) { // [max(kp,kp-q),+inf)
          double temp = 0., r = k0-p-q;
          double v=k0-q;
          temp += lga( qp*qm/(kp*km) )*F_123(q,p,r)/v ;   // 4C
          return temp;                                    // 4B = 0
  },  max(kp,kp-q)  )(y); },  km,kp    )(x); //*/
  _4B +=
  remap([&](double q, double qd) { // [km,0]
      qm = km-q;
      return remap([&](double p, double pd) { // [kp,kp-q]
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          pp = (fabs(pp)<1e-1*fabs(q)) ? pd : pp;
          double l=k0-p, v=k0-q;
          temp += lga( km*q/(pp*r) )*F_123(p,q,r)/l ;     // 4B'
          temp += lga( qm*p/(kp*r) )*F_123(q,p,r)/v ;     // 4C'
          return temp;
  },  kp,kp-q  )(y); },  km,0.    )(x); //*/\
  return _4B;
  _4D =
  remap([&](double q, double qd) { // [0,kp]
      return remap([&](double p, double pd) { // [kp-q,kp]
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          pp = (fabs(pp)<1e-1*fabs(q)) ? pd : pp;
          double l=k0-p;
          temp += lga( pp*r/(q*km) )*F_123(p,q,r)/l ;     // 4D
          return temp;
  },  kp-q,kp  )(y); },  0. ,kp    )(x); //*/\
  return _4D;

  //  #] from '4'\
      #[ from '5'

  _5A =
  remap([&](double q, double qd) { // (-inf,km]
      return remap([&](double p, double pd) { // (-inf,km]
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          double l=k0-p;
          temp += lga( q*pp/(r*km) )*F_123(p,q,r)/l ;     // 5A
          return temp;
  },  km  )(y); },  km    )(x); //*/\
  return _5A;
  _5B =
  remap([&](double q, double qd) { // [km,0]
      return remap([&](double p, double pd) { // [km,km-q]
          double temp = 0., r = k0-p-q;
          pm = km-p;
          pm = (fabs(pm)<1e-1*fabs(q)) ? pd : pm;
          double l=k0-p;
          temp += lga( pm*r/(kp*q) )*F_123(p,q,r)/l ;     // 5B
          return temp;
  },  km,km-q  )(y); },  km,0.    )(x); //*/\
  return _5B;

  //  #] from '5'\
      #[ from '6'

  _6A =
  remap([&](double p, double pd) { // [k,+inf)
      return remap([&](double q, double qd) { // [kp-p,km]
          double temp = 0., r = k0-p-q;
          pp = kp-p;
          qp = kp-q;
          double l=k0-p, v=k0-q;
          temp += lga( q*pp/(km*r) )*F_123(p,q,r)/l ;     // 6A
          temp += lga( p*qp/(km*r) )*F_123(q,p,r)/v ;     // 7A
          return temp;
  },  kp-p,km  )(y); },  k    )(x); //*/\
  return _6A;
  _6B =
  remap([&](double pq, double pqd) { // [k-km,+inf)
      return remap([&](double r, double rd) { // [km,-km]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          double l=k0-p, v=k0-q;
          temp += lga( q/r )*F_123(p,q,r)/l ;             // 6B
          temp += lga( p/r )*F_123(q,p,r)/v ;             // 7B
          return temp;
  },  km,-km  )(y); },  k-km    )(x); //*/
  _6B +=
  remap([&](double pq, double pqd) { // [k-pq,pq-k]
      return remap([&](double r, double rd) { // [k,k-km]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          double l=k0-p, v=k0-q;
          temp += lga( q/r )*F_123(p,q,r)/l ;             // 6B
          temp += lga( p/r )*F_123(q,p,r)/v ;             // 7B
          return temp;
  },  k-pq,pq-k  )(y); },  k,k-km    )(x); //*/
  _6B +=
  remap([&](double p, double pd) { // [kp,+inf)
      return remap([&](double q, double qd) { // [km-p,km+k0-p]
          double temp = 0., r = k0-p-q;
          double l=k0-p, v=k0-q;
          temp += lga( q/r )*F_123(p,q,r)/l ;             // 6B
          temp += lga( p/r )*F_123(q,p,r)/v ;             // 7B
          return 2.*temp;
  },  km-p, km+k0-p  )(y); },  kp    )(x); //*/\
  return _6B;
  _6C =
  remap([&](double q, double qd) { // (-inf,-k]
      return remap([&](double p, double pd) { // [kp,km-q]
          double temp = 0., r = k0-p-q;
          pm = km-p;
          qm = km-q;
          double l=k0-p, v=k0-q;
          temp += lga( q*pm/(kp*r) )*F_123(p,q,r)/l ; // 6C
          temp += lga( p*qm/(kp*r) )*F_123(q,p,r)/v ; // 7C
          return temp;
  },  kp,km-q  )(y); },  -k    )(x); //*/\
  return _6C;

  //  #] from '6',\
      #[ from virtual corrections

  _Virt =
  remap([&](double p, double pd) { // [kp,+inf)
      double tempp = 0.;
      double l=k0-p;
      pp = kp-p; pm = km-p;
      pp = (fabs(pp)<1e-1*k0) ? pd : pp; // q=kp (!)
      tempp += remap([&](double q, double qd) { // [0,2l]
          double temp = 0., r=k0-p-q, rp=k0-p+q;

          temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*pow(q/k0,n) ;
          temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*pow(-q/k0,n) ;
          return temp;
      },  0.,2.*fabs(k0-p) )(y); 
      tempp += remap([&](double q, double qd) { // [2l,+inf)
          double temp = 0., r=k0-p-q, rp=k0-p+q;

          temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*pow(q/k0,n) ;
          temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*pow(-q/k0,n) ;
          return temp;
      },  2.*fabs(k0-p) )(y); 

      return F_14(p,l)*(
                      ( lga(4.*pm*pp/(k*k))+lga(l*l/K2)
                    + ( (m+n==1) ? 3.5 : 3. ) // 3.5 if m=1 !!!!
                    )*pow( .5*l/k0, n )
                    + 2.*s3*tempp/l
                    );
  }, kp    )(x); //*/
  _Virt +=
  remap([&](double p, double pd) { // (-inf,km]
      double tempp = 0.;
      double l=k0-p;
      pp = kp-p; pm = km-p;
      pm = (fabs(pm)<1e-1*k0) ? pd : pm; // q=km (!)
      tempp += remap([&](double q, double qd) { // [0,2l]
          double temp = 0., r=k0-p-q, rp=k0-p+q;

          temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*pow(q/k0,n) ;
          temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*pow(-q/k0,n) ;
          return temp;
      },  0.,2.*fabs(k0-p) )(y); 
      tempp += remap([&](double q, double qd) { // [2l,+inf)
          double temp = 0., r=k0-p-q, rp=k0-p+q;

          temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*pow(q/k0,n) ;
          temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*pow(-q/k0,n) ;
          return temp;
      },  2.*fabs(k0-p) )(y); 

      return F_14(p,l)*(
                      ( lga(4.*pm*pp/(k*k))+lga(l*l/K2)
                    + ( (m+n==1) ? 3.5 : 3. ) // 3.5 if m=1 !!!!
                    )*pow( .5*l/k0, n )
                    + 2.*s3*tempp/l
                    );
  }, km    )(x); //*/\

  res += _2A + _2B + _4A + _4B + _4D + _5A + _5B + _6A + _6B + _6C + _Virt;

  }
    // still don't have a good way to cater for NaNs
  if ( isinf(res)||isnan(res) ) { return 0.;}
  else { return .25*fabs(K2)*res/(k); }
}
