#include "core.hh"
#include "quad.hh"
#include "map.hh"

/*--------------------------------------------------------------------*/

struct rhoStar : Master {
  double integrand(double,double); // supported on [0,1]x[0,1]

  double F_123(double,double,double); // Cut: f1f2f3/f0   }- real
  double F_14(double,double);         //        f1f4/f0   }- virtual

  double g_(double,double);            // non-log function

  double eval();
  rhoStar(int _s[3]) : Master(0,0,_s) { type=7; }
};
// function for MAIN
Master* _Star(int m, int n, int s[3]) {
  Master *R =  new rhoStar(s); return R;
}

/*--------------------------------------------------------------------*/
// all the thermal weights for cuts from this topology

double rhoStar::F_123(double p, double q, double r) { 
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

double rhoStar::F_14(double p, double l) {
  double res;

  int s1 = (this->s)[1],
      s4 = (this->s)[4];

  int _m = this->m;

  double fp=f(p,s1), fl=f(l,s4);

  res = ( ((double) s1*s4)*exp(k0)-1. )*fp*fl;
  res*= pow(k0,-_m)*pow(p,_m)*sgn(km);

  return res;
}

double rhoStar::g_(double p, double q) {
  double r=k0-p-q, l=k0-p, 
         res;
  double A, B; // min & max resp.
  A = (2.*q*k0 - K2)/(2.*q*k); 
  B = (2.*q*k0 - K2)/(2.*q*k); 
  double tmp =2.*p*r/(k*q);
  (tmp<0) ? A+=tmp : B+=tmp; // ensure A < B
  A = max(-1.,A);
  B = min( 1.,B);

  double gu = (k0*k0+k*k-2.*k0*p-2.*k*l*B); gu*=gu;
  double gl = (k0*k0+k*k-2.*k0*p-2.*k*l*A); gl*=gl;

  res = sqrt(gu)-sqrt(gl);

  res *= -sgn(p);
  if (A> 1.) { res=0.; }    // just a failsafe, not nec
  if (B<-1.) { res=0.; }
  if (A> B ) { res=0.; }

  bool reg=((km>p)||(kp<p));
  if (( (!reg)&&(k0>k) ) || ( (reg)&&(k0<k) )) {
    res += 2.*sgn(km)*( K2 - 2.*l*k0);
  }

  return res*q;
}

/*--------------------------------------------------------------------*/
// evaluation step, including the OPE

double rhoStar::eval() {
  double res, err;

  double a1=I(0,(this->s)[1]), b1=I(2,(this->s)[1]), // tadpole ints
         a2=I(0,(this->s)[2]), b2=I(2,(this->s)[2]),
         a3=I(0,(this->s)[3]), b3=I(2,(this->s)[3]),
         a4=I(0,(this->s)[4]), b4=I(2,(this->s)[4]);

  if ( m==0 && n==0 ) { // K.Q
    (this->OPE).T0 = -K2*5./4.;
    (this->OPE).T2 = +( a1 )*.25*OOFP;
    (this->OPE).T4 = +( 3.*(b1-b2)+b3 )*(k0*k0+k*k/3.)/SQR(K2)/6.*OOFP;
  } 
  //return (this->OPE)();

  double epsabs = 1e-4, epsrel = 0;
  //if (k0>20.) { epsabs*=.1; epsrel*=.1; }
  size_t limit = 1e5;
  gsl_set_error_handler_off(); // live on the edge.

  quad wsp1(limit);
  quad wsp2(limit);

  auto outer = make_gsl_function( [&](double x) {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
        return (this->integrand)(x,y);
        } );
    gsl_integration_qag(inner, .0+1e-16,1., epsabs, epsrel, 
                       limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag(  outer, .0+1e-16,1., epsabs, 2.*epsrel, 
                        limit, 6, wsp2, &res, &err  );

  return (( res*(kp/km) ))*CUBE(OOFP); // mult. by (km/kp) in integrand
                                       // to soften LC divergence
}

/*--------------------------------------------------------------------*/

double rhoStar::integrand(double x, double y) 
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
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      qm = km-q; qp = kp-q;
      if (k0>3.*k) {
        temp += ( g_(p,q) + q*K2*lga( kp*pp/(km*pm) ) )*F_123(p,q,r)/l2 ;
        temp += ( g_(q,p) + p*K2*lga( kp*qp/(km*qm) ) )*F_123(q,p,r)/v2 ; // 1A
      } else {
        temp += ( q*K2*lga( pp*r/(km*q) ) )*F_123(p,q,r)/l2 ;
        temp += ( p*K2*lga( qp*r/(km*p) ) )*F_123(q,p,r)/v2 ; // 1A'
        temp *= -1.; // (signed integration region)
        // p -> -p and q -> k0-q (missing part of 2D, 3D)
        temp += ( g_(l,-q) )*F_123(l,-q,p+q)/SQR(p) ;
        temp += ( g_(v,-p) )*F_123(v,-p,p+q)/SQR(q) ; // 
        temp += ( g_(-p,v) + v*K2*lga( km*kp/((kp+p)*(km+p)) ) )*
                F_123(-p,v,p+q)/SQR(k0+p) ;
        temp += ( g_(-q,l) + l*K2*lga( km*kp/((kp+q)*(km+q)) ) )*
                F_123(-q,l,p+q)/SQR(k0+q) ;
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
      qm = km-q;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( kp*q/(r*pm) ) )*F_123(p,q,r)/l2 ;
      temp += ( g_(q,p) + p*K2*lga( kp*p/(r*qm) ) )*F_123(q,p,r)/v2 ; // 1B
      return temp;
  },  km-p,min(p,kp-p) )(y); }, .5*km,km   )(x);//*/
  //\
  return _1B;
  //
  _1D = 
  remap([&](double p, double pd) { // [max(k,km),kp]
    return remap([&](double q, double qd) { // [kp-p,km]
      qm = km-q;
      pp = kp-p;
      //\
      Padding for log singularities:
      double h=fabs(km-kp+p);
      qm = (fabs(qm)<1e-1*h) ? qd : qm;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( pp*r/(km*q) ) )*F_123(p,q,r)/l2 ; // 1D
      temp += ( g_(q,p) + p*K2*lga( kp*r/(p*qm) ) )*F_123(q,p,r)/v2 ; // 1C
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
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      qm = km-q; qp = kp-q;
      temp += ( g_(p,q) + q*K2*lga( km*pm/(kp*pp) ) )*F_123(p,q,r)/l2 ; // 2A
      temp += ( g_(q,p) + p*K2*lga( km*qm/(kp*qp) ) )*F_123(q,p,r)/v2 ; // 3A
      return temp;
  },  km-p  )(y); },  kp,k0+km )(x); //*/
  _2A +=
  remap([&](double p, double pd) { // [k0+km,+inf)
    pm = km-p; pp = kp-p;
    return remap([&](double q, double qd) { // (-inf,km-p]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      qm = km-q; qp = kp-q;
      temp += ( g_(p,q) + q*K2*lga( km*pm/(kp*pp) ) )*F_123(p,q,r)/l2 ; // 2A
      temp += ( g_(q,p) + p*K2*lga( km*qm/(kp*qp) ) )*F_123(q,p,r)/v2 ; // 3A
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
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      qp = kp-q;
      temp += ( g_(p,q) + q*K2*lga( km*q/(pp*r) ) )*F_123(p,q,r)/l2 ; // 2B
      temp += ( g_(q,p) + p*K2*lga( km*p/(qp*r) ) )*F_123(q,p,r)/v2 ; // 3B
      return temp;
  },  km-p,kp-p  )(y); },  kp,k0+km )(x); //*/
  _2B +=
  remap([&](double p, double pd) { // [k0+km,+inf)
    pp = kp-p;
    return remap([&](double q, double qd) { // [km-p,kp-p]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      qp = kp-q;
      temp += ( g_(p,q) + q*K2*lga( km*q/(pp*r) ) )*F_123(p,q,r)/l2 ; // 2B
      temp += ( g_(q,p) + p*K2*lga( km*p/(qp*r) ) )*F_123(q,p,r)/v2 ; // 3B
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
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      qp = kp-q;
      temp += ( g_(p,q) + q*K2*lga( pm*r/(kp*q) ) )*F_123(p,q,r)/l2 ; // 2B
      temp += ( g_(q,p) + p*K2*lga( km*r/(qp*p) ) )*F_123(q,p,r)/v2 ; // 3B
      return temp;
  },  km-p  )(y); },  km,kp )(x); //*/
  //\
  return _2C;
  _2D =
  remap([&](double q, double qd) { // [-min(k,km),0]
    qm = km-q; qp = kp-q;
    return remap([&](double p, double pd) { // [km-q,kp]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q)  )*F_123(p,q,r)/l2 ;// 2D
      temp += ( g_(l,-q) )*F_123(l,-q,p+q)/SQR(p) ; // 1C'
      temp += ( g_(q,p) + p*K2*lga( km*kp/(qp*qm) ) )*F_123(q,p,r)/v2 ;// 3D \
                                                               (+part from 1A)
      temp += ( g_(-q,l) + l*K2*lga( km*kp/((kp+q)*(km+q)) ) )
                                           *F_123(-q,l,p+q)/SQR(k0+q) ; // 1D'
      return temp;
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
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( km*pm/(kp*pp) ) )*F_123(p,q,r)/l2 ;//    \
                                                                 no reflection
      return temp;
  },  kp  )(y); }, kp,k0+km    )(x); //*/
  _4A +=
  remap([&](double p, double pd) { // [k0+km,+inf)
    pm = km-p; pp = kp-p;
    return remap([&](double q, double qd) { // [kp,+inf)
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( km*pm/(kp*pp) ) )*F_123(p,q,r)/l2 ;//    \
                                                                 no reflection
      return temp;
  },  kp  )(y); },  k0+km    )(x); //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) { // [km,kp]
    qm = km-q;
    qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=kp (!)
    return remap([&](double p, double pd) { // [kp,k0+km]
      pp = kp-p;
      // the sing. @ p=k0 is at the midpoint of the interval
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( km*r/(q*pp) ) )*F_123(p,q,r)/l2 ;// 4B
      temp += ( g_(q,p) + p*K2*lga( qm*r/(p*kp) ) )*F_123(q,p,r)/v2 ;// 4C
      return temp;
  },  kp,k0+km  )(y); }, km,kp    )(x); //*/
  _4B +=
  remap([&](double q, double qd) { // [km,kp]
    qm = km-q;
    qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=kp (!)
    return remap([&](double p, double pd) { // [kp,k0+km]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      pp = kp-p;
      temp += ( g_(p,q) + q*K2*lga( km*r/(q*pp) ) )*F_123(p,q,r)/l2 ;// 4B
      temp += ( g_(q,p) + p*K2*lga( qm*r/(p*kp) ) )*F_123(q,p,r)/v2 ;// 4C
      return temp;
  },  k0+km  )(y); }, km,kp    )(x); //*/
  //\
  return _4B;
  _4D =
  remap([&](double pq, double pqd) { // [0,k]
    return remap([&](double r, double rd) { // [-k+p-q,k-p+q]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( .5*g_(p,q) + q*K2*lga( r/q ) )*F_123(p,q,r)/l2 ;// 4D+1B+1B'
      temp += ( .5*g_(q,p) + p*K2*lga( r/p ) )*F_123(q,p,r)/v2 ;// (reflected)
      return temp; // no .5, from log(...^2)
  },  -k+pq,k-pq )(y); }, 0.,k    )(x); //*/
  //\
  return _4D;

  //  #] from '4',\
      #[ from virtual corrections

  _Virt =
  remap([&](double p, double pd) { // [km,kp]
      double tempp = 0.;
      double l=k0-p, l2=l*l;
      pp = kp-p; pm = km-p;
      pp = (fabs(pp)<1e-1*k) ? pd : pp; // q=kp (!)
      pm = (fabs(pm)<1e-1*k) ? pd : pm; // q=km (!)
      tempp += remap([&](double q, double qd) { // [0,2l]
        double temp = 0., r=k0-p-q, rp=k0-p+q;

        temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*( +q) ;
        temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*( -q )  ;
        return temp;
      },  0.,2.*(k0-p) )(y); 
      tempp += remap([&](double q, double qd) { // [2l,+inf)
        double temp = 0., r=k0-p-q, rp=k0-p+q;

        temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*( +q) ;
        temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*( -q) ;
        return temp;
      },  2.*(k0-p) )(y); 

      return K2*F_14(p,l)*(
                      .5*( lga(4.*pm*pp/(k*k))+lga(l2/K2) )
                    + 7./4.
                    - (K2 - 2.*k0*l)/K2
                    + 2.*tempp/l2
                    );
  }, km,kp    )(x); //*/
  //\
  return .25*_Virt/k;

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
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( pp*pm/(kp*km) ) )*F_123(p,q,r)/l2 ;// 2A
      temp += ( g_(q,p) )*F_123(q,p,r)/v2 ;// 3A
      return temp;
  },  km  )(y); },  km,0.    )(x);
  _2A +=
  remap([&](double p, double pd) { // [0,k0+km]
    pm = km-p, pp = kp-p;
    return remap([&](double q, double qd) { // (-inf,km-p]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( pp*pm/(kp*km) ) )*F_123(p,q,r)/l2 ;// 2A
      temp += ( g_(q,p) )*F_123(q,p,r)/v2 ;// 3A
      return temp;
  }, km-p  )(y); },  0.,k0+km    )(x);
  _2A +=
  remap([&](double p, double pd) { // [k0+km,kp]
    pm = km-p, pp = kp-p;
    return remap([&](double q, double qd) { // (-inf,km-p]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( pp*pm/(kp*km) ) )*F_123(p,q,r)/l2 ;// 2A
      temp += ( g_(q,p) )*F_123(q,p,r)/v2 ;// 3A
      return temp;
  },  km-p  )(y); },  k0+km,kp    )(x); //*/
  //\
  return _2A;
  _2B  =
  remap([&](double p, double pd) { // [0,k0+km]
    pp = kp-p;
    return remap([&](double q, double qd) { // [km-p,km]
      qm = km-q;
      qm = (fabs(qm)<1e-1*p) ? qd : qm;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( q*pp/(r*km) ) )*F_123(p,q,r)/l2 ;// 2B
      temp += ( g_(q,p) + p*K2*lga( p*kp/(r*qm) ) )*F_123(q,p,r)/v2 ;// 3B
      return temp;
  },  km-p,km  )(y); },  0.,k0+km    )(x);
  _2B +=
  remap([&](double p, double pd) { // [k0+km,kp]
    pp = kp-p;
    return remap([&](double q, double qd) { // [km-p,km]
      qm = km-q;
      qm = (fabs(qm)<1e-1*p) ? qd : qm;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( q*pp/(r*km) ) )*F_123(p,q,r)/l2 ;// 2B
      temp += ( g_(q,p) + p*K2*lga( p*kp/(r*qm) ) )*F_123(q,p,r)/v2 ;// 3B
      return temp;
  },   km-p,km  )(y); },  k0+km,kp    )(x); //*/
  //\
  return _2B;

  //  #] from '2'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) { // [kp,+inf)
    pm = km-p;
    return remap([&](double q, double qd) { // [kp,+inf)
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l;
      temp += ( g_(p,q) + q*K2*lga( q*pm/(r*kp) ) )*F_123(p,q,r)/l2 ;// 4A
      return temp;
  },  kp  )(y); },  kp    )(x); //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) { // [km,kp]
    qp = kp-q; qm = km-q;
    return remap([&](double p, double pd) { // [max(kp,kp-q),+inf)
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) )*F_123(p,q,r)/l2 ;// 4B
      temp += ( g_(q,p) + p*K2*lga( qp*qm/(kp*km) ) )*F_123(q,p,r)/v2 ;// 4C
      return temp;
  },  max(kp,kp-q)  )(y); },  km,kp    )(x); //*/
  _4B +=
  remap([&](double q, double qd) { // [km,0]
    qm = km-q;
    return remap([&](double p, double pd) { // [kp,kp-q]
      pp = kp-p;
      pp = (fabs(pp)<1e-1*fabs(q)) ? pd : pp;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( km*q/(pp*r) ) )*F_123(p,q,r)/l2 ;// 4B'
      temp += ( g_(q,p) + p*K2*lga( qm*p/(kp*r) ) )*F_123(q,p,r)/v2 ;// 4C'
      return temp;
  },  kp,kp-q  )(y); },  km,0.    )(x); //*/
  //\
  return _4B;
  _4D =
  remap([&](double q, double qd) { // [0,kp]
    return remap([&](double p, double pd) { // [kp-q,kp]
      pp = kp-p;
      pp = (fabs(pp)<1e-1*fabs(q)) ? pd : pp;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( pp*r/(q*km) ) )*F_123(p,q,r)/l2 ;// 4D
      return temp;
  },  kp-q,kp  )(y); },  0. ,kp    )(x); //*/
  //\
  return _4D;

  //  #] from '4'\
      #[ from '5'

  _5A =
  remap([&](double q, double qd) { // (-inf,km]
    return remap([&](double p, double pd) { // (-inf,km]
      pp = kp-p;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( q*pp/(r*km) ) )*F_123(p,q,r)/l2 ;// 5A
      return temp;
  },  km  )(y); },  km    )(x); //*/
  //\
  return _5A;
  _5B =
  remap([&](double q, double qd) { // [km,0]
    return remap([&](double p, double pd) { // [km,km-q]
      pm = km-p;
      pm = (fabs(pm)<1e-1*fabs(q)) ? pd : pm;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( pm*r/(kp*q) ) )*F_123(p,q,r)/l2 ;// 5B
      return temp;
  },  km,km-q  )(y); },  km,0.    )(x); //*/
  //\
  return _5B;

  //  #] from '5'\
      #[ from '6'

  _6A =
  remap([&](double p, double pd) { // [k,+inf)
    return remap([&](double q, double qd) { // [kp-p,km]
      pp = kp-p;
      qp = kp-q;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( q*pp/(km*r) ) )*F_123(p,q,r)/l2 ;// 6A
      temp += ( g_(q,p) + p*K2*lga( p*qp/(km*r) ) )*F_123(q,p,r)/v2 ;// 7A
      return temp;
  },  kp-p,km  )(y); },  k    )(x); //*/
  //\
  return _6A;
  _6B =
  remap([&](double pq, double pqd) { // [k-km,+inf)
    return remap([&](double r, double rd) { // [km,-km]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( .5*g_(p,q) + q*K2*lga( q/r ) )*F_123(p,q,r)/l2 ;// 6B
      temp += ( .5*g_(q,p) + p*K2*lga( p/r ) )*F_123(q,p,r)/v2 ;// 7B
      return temp;
  },  km,-km  )(y); },  k-km    )(x); //*/
  _6B +=
  remap([&](double pq, double pqd) { // [k-pq,pq-k]
    return remap([&](double r, double rd) { // [k,k-km]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( .5*g_(p,q) + q*K2*lga( q/r ) )*F_123(p,q,r)/l2 ;// 6B
      temp += ( .5*g_(q,p) + p*K2*lga( p/r ) )*F_123(q,p,r)/v2 ;// 7B
      return temp;
  },  k-pq,pq-k  )(y); },  k,k-km    )(x); //*/
  _6B +=
  remap([&](double p, double pd) { // [kp,+inf)
    return remap([&](double q, double qd) { // [km-p,km+k0-p]
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( .5*g_(p,q) + q*K2*lga( q/r ) )*F_123(p,q,r)/l2 ;// 6B
      temp += ( .5*g_(q,p) + p*K2*lga( p/r ) )*F_123(q,p,r)/v2 ;// 7B
      return 2.*temp;
  },  km-p, km+k0-p  )(y); },  kp    )(x); //*/\
  return _6B;
  _6C =
  remap([&](double q, double qd) { // (-inf,-k]
    return remap([&](double p, double pd) { // [kp,km-q]
      pm = km-p;
      qm = km-q;
      double temp = 0., r = k0-p-q;
      double l=k0-p, l2=l*l, v=k0-q, v2=v*v;
      temp += ( g_(p,q) + q*K2*lga( q*pm/(kp*r) ) )*F_123(p,q,r)/l2 ;// 6C
      temp += ( g_(q,p) + p*K2*lga( p*qm/(kp*r) ) )*F_123(q,p,r)/v2 ;// 7C
      return temp;
  },  kp,km-q  )(y); },  -k    )(x); //*/\
  return _6C;

  //  #] from '6',\
      #[ from virtual corrections

  _Virt =
  remap([&](double p, double pd) { // [kp,+inf)
    double tempp = 0.;
    double l=k0-p, l2=l*l;
    pp = kp-p; pm = km-p;
    pp = (fabs(pp)<1e-1*k0) ? pd : pp; // q=kp (!)
    tempp += remap([&](double q, double qd) { // [0,2l]
      double temp = 0., r=k0-p-q, rp=k0-p+q;
      temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*( +q ) ;
      temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*( -q ) ;
      return temp;
    },  0.,2.*fabs(k0-p) )(y); 
    tempp += remap([&](double q, double qd) { // [2l,+inf)
      double temp = 0., r=k0-p-q, rp=k0-p+q;
      temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*( +q ) ;
      temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*( -q ) ;
      return temp;
    },  2.*fabs(k0-p) )(y);
    //
    return K2*F_14(p,l)*(
                          + .5*( lga(4.*pm*pp/(k*k))+lga(l2/K2) )
                          + 7./4.
                          - (K2 - 2.*k0*l)/K2
                          + 2.*tempp/l2
                        );
  }, kp    )(x); //*/
  _Virt +=
  remap([&](double p, double pd) { // (-inf,km]
    double tempp = 0.;
    double l=k0-p, l2=l*l;
    pp = kp-p; pm = km-p;
    pm = (fabs(pm)<1e-1*k0) ? pd : pm; // q=km (!)
    tempp += remap([&](double q, double qd) { // [0,2l]
        double temp = 0., r=k0-p-q, rp=k0-p+q;
        temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*( +q ) ;
        temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*( -q ) ;
        return temp;
    },  0.,2.*fabs(k0-p) )(y); 
    tempp += remap([&](double q, double qd) { // [2l,+inf)
        double temp = 0., r=k0-p-q, rp=k0-p+q;
        temp += lga( q/r )*(sgn(r)*f(fabs(r),s3))*( +q ) ;
        temp += lga( q/rp)*(sgn(rp)*f(fabs(rp),s3))*( -q ) ;
        return temp;
    },  2.*fabs(k0-p) )(y); 
    //
    return K2*F_14(p,l)*(
                          + .5*( lga(4.*pm*pp/(k*k))+lga(l2/K2) )
                          + 7./4.
                          - (K2 - 2.*k0*l)/K2
                          + 2.*tempp/l2
                        );
  }, km    )(x); //*/
  //\
  return _Virt;

  res +=
  _2A + _2B + _4A + _4B + _4D + _5A + _5B + _6A + _6B + _6C + _Virt;

  }
    // still don't have a good way to cater for NaNs
  if ( isinf(res)||isnan(res) ) { return 0.;}
  else { return .25*(km/kp)*res/(k); }
}

