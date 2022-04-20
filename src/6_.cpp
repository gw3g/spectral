#include "core.hh"
#include "quad.hh"
#include "map.hh"

double checker(double x) {
  if ( isinf(x)||isnan(x) ) { return 0.;}
  else return x;
}

using namespace std;
double E = 1e2; // control param for large momenta: used to \
                   switch between expanded or full integrand\
                   See: L_div(...) and L_int(...)

/*--------------------------------------------------------------------*/

struct rho11111 : Master {
  double integrand(double,double); // supported on [0,1]x[0,1]

  double F_123(double,double,double); // Cut: f1f2f3/f0   \_ real
  double F_345(double,double,double); //      f3f4f5/f0   /
  double F_14(double,double);         //        f1f4/f0   \_ virtual
  double F_25(double,double);         //        f2f5/f0   /

  double eval();

  rho11111(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=6; }
};
// function for MAIN
Master* _11111(int m, int n, int s[3]) {
  Master *R =  new rho11111(m,n,s); return R;
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

  double fp=f(p-MOT1,s1), fq=f(q-MOT2,s2), fr= f(r-MOT3,s3);
  //double fp = f(p,s1) ? p>0 : -1.-f(-p,s1);
  //double fq = f(q,s2) ? q>0 : -1.-f(-q,s2);
  //double fr = f(r,s3) ? r>0 : -1.-f(-r,s3);

  res = ( ((double) s1*s2*s3)*exp(k0)-1. )*fp*fq*fr;
  //res = 1. + fp + fq + fr + fp*fq + fp*fr + fq*fr;
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

  double fr=f(r-MOT3,s3), fl=f(l-MOT4,s4), fv= f(v-MOT5,s5);
  //double fr = f(r,s3) ? r>0 : -1.-f(-r,s3);
  //double fl = f(l,s4) ? l>0 : -1.-f(-l,s4);
  //double fv = f(v,s5) ? v>0 : -1.-f(-v,s5);

  res = ( ((double) s3*s4*s5)*exp(k0)-1. )*fr*fl*fv;
  //res = 1. + fl + fv + fr + fl*fv + fl*fr + fv*fr;
  res*= pow(k0,-_m-_n)*pow(k0-l,_m)*pow(k0-v,_n) ;

  return res;
}

double rho11111::F_14(double p, double l) {
  double res;

  int s1 = (this->s)[1],
      s4 = (this->s)[4];

  int _m = this->m;

  double fp=f(p-MOT1,s1), fl=f(l-MOT4,s4);
  //double fp = f(p,s1) ? p>0 : -1.-f(-p,s1);
  //double fl = f(l,s4) ? l>0 : -1.-f(-l,s4);

  res = ( ((double) s1*s4)*exp(k0)-1. )*fp*fl;
  //res = 1. + fp + fl;
  res*= pow(k0,-_m)*pow(p,_m)*sgn(km);

  return res;
}

double rho11111::F_25(double q, double v) {
  double res;

  int s2 = (this->s)[2],
      s5 = (this->s)[5];

  int _n = this->n;

  double fq=f(q-MOT2,s2), fv=f(v-MOT5,s5);
  //double fq = f(q,s2) ? q>0 : -1.-f(-q,s2);
  //double fv = f(v,s5) ? v>0 : -1.-f(-v,s5);

  res = ( ((double) s2*s5)*exp(k0)-1. )*fq*fv;
  //res = 1. + fq + fv;
  res*= pow(k0,-_n)*pow(q,_n)*sgn(km);

  return res;
}

/*--------------------------------------------------------------------*/
// evaluation step, including the OPE

double rho11111::eval() 
{
  double res, err;

  double a1=I(0,(this->s)[1],MOT1), b1=I(2,(this->s)[1],MOT1), // tadpole ints
         a2=I(0,(this->s)[2],MOT2), b2=I(2,(this->s)[2],MOT2),
         a3=I(0,(this->s)[3],MOT3), b3=I(2,(this->s)[3],MOT3),
         a4=I(0,(this->s)[4],MOT4), b4=I(2,(this->s)[4],MOT4),
         a5=I(0,(this->s)[5],MOT5), b5=I(2,(this->s)[5],MOT5);

  if ( m==0 && n==0 ) { // (0)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 = -( a1+a2+2.*a3+a4+a5 )*.25*OOFP;
    (this->OPE).T4 = -( (b1+b2+b4+b5)*11.+b3*6.)*(k0*k0+k*k/3.)/SQR(K2)/6.*OOFP;
  } else
  if ( m==1 && n==0 ) { // (1)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 = -( a2+a3+a4 )*.25*OOFP*k0;
    (this->OPE).T4 = -( ( (b2+b4)*11.+b3*3.)*(k0*k0+k*k/3.) +
                        ( (b1-b4)*9.-(b2-b5)*5. )*K2*.5     )*k0/SQR(K2)/6.*OOFP;
  } else
  if ( m==0 && n==1 ) { // (0,1)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 = -( a1+a3+a5 )*.25*OOFP*k0;
    (this->OPE).T4 = -( ( (b1+b5)*11.+b3*3.)*(k0*k0+k*k/3.) +
                        ( (b2-b5)*9.-(b1-b4)*5. )*K2*.5     )*k0/SQR(K2)/6.*OOFP;
  } else
  if ( m==1 && n==1 ) { // (1,1)
    (this->OPE).T0 = 1./16.;
    (this->OPE).T2 = 0.;
    (this->OPE).T4 = +(( (b1+b2+b4+b5)*3.-b3*2.      )*.5*K2 -
                       ( (b1+b2)*9.+b3*2.+(b4+b5)*5. )*k0*k0 )/K2/12.*OOFP;
  } else
  if ( m==2 && n==0 ) { // (2)
    (this->OPE).T0 = -.5*( +11./16. -.5 );
    (this->OPE).T2 = -( (a2+a3+a4)*k0*k0-(a2+a5)*K2*.25 )*.25*OOFP;
    (this->OPE).T4 = -(((b1+b4)*3.+(b2+b3+b5)*2.)*.5*K2-
                       (b2*7.+b3+b4*9.+b5*2.)*k0*k0+
                       ((b2+b4)*11.+3.*b3 )*k0*k0*(k0*k0+k*k/3.)/K2  )/K2/6.*OOFP;
  } else {
    cerr << "Case: (m,n)=("<< m << ","<<n<<") out of bounds!\n";
    return 0.;
  }

  // consistent dimensions
  (this->OPE).T2 /= SQR(K2);
  (this->OPE).T4 /= SQR(K2);
  //return (this->OPE)(); 

  // Quadrature step! --
  double epsabs = 1e-4, epsrel = 0;
  size_t limit = 5e2;

  gsl_set_error_handler_off();
  quad wsp1(limit);
  quad wsp2(limit);

  if ( m==2 && n==0 ) { epsabs = 2e-3; } // lower tol for m=2 master

  auto outer = make_gsl_function( [&](double x) 
  {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
          return (this->integrand)(x,y);
        } );
    gsl_integration_qag( inner, .0+1e-13,1., epsabs, epsrel,
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag( outer, .0+1e-13,1., epsabs, epsrel*2,
                       limit, 6, wsp2, &res, &err  );//*/

  // Rescalings --
  if ( m==2 && n==0 ) {
    res = res*(kp/km)*CUBE(OOFP)/(K2); return res;
  }
  if ( m==1 && n==1 ) {
  res = res*(kp/km)*CUBE(OOFP)/K2; return res;
  }
  res = res*(kp/km)*pow(k0,m+n)*CUBE(OOFP)/SQR(K2); return res;
}

/*--------------------------------------------------------------------*/
// integrand manipulations

inline double L_div(double p, double l, double pp, double pm) {
  //double pp = kp-p, pm = km-p, l=k0-p;
  double res = lga( pp*pm/(l*p) );//-K2/(4.*p*p);
  if (fabs(p/kp)>E) {
    res = 0.; double Rp=SQR(kp/p),
                     Rm=SQR(km/p),
                     R0=SQR(k0/p);

    for (int n=3;n<100;n++) { // series expansion, p -> +inf
      Rp *= kp/p;
      Rm *= km/p;
      R0 *= k0/p;
      res -= ( Rp+Rm-R0 )/((double)n);
    }; // --> first term omitted
       //     ( it is K2/(4*p^2) )

  }
  return res;
}

inline double L_int(int j, double a, double p) {
  // indefinite integrals of the form
  //
  //  ∫dp (p/k0)^j (p^2*(a-p))^{-1} ; _/ just indefinite result,
  //                                   \ with p=+inf subtracted off]
  double temp = lga(p/(p-a)), res=0;
  if (j==0) {
    res = ( temp/a - 1./p )/a;
  } else
  if (j==1) {
    res = temp/(a*k0);
  } else
  if (j==2) { //return 0.;
    res =  ( - 5.5
             + lga(.25*k*k/((km-a)*(kp-a)))
             + 2. // finite part, consistent w/ cut-off
             - lga( SQR(p-a)/K2 )                 )/(2.*k0*k0);
    // div = 1/eps + 2*log(mu2/K2) + 11/2.
  }
  return res*(K2/16.); // Why 16? Long story ...
}

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
  remap([&](double p, double pd) {                      // p=[.5kp,km]
    return remap([&](double q, double qd) {             // q=[kp-p,p]
      double temp = 0., r=k0-p-q;
      pm = km-p; pp = kp-p;
      qm = km-q; qp = kp-q;
      pm = (fabs(pm)<1e-1*fabs(km-k)) ? pd : pm;        // p=km (!)
      if (k0>3.*k) {
        temp +=  lga(pp/pm)*(                           // refl. log
                              F_123(p,q,r)              // 1A[p,q]
                            + F_123(q,p,r)              // 1A[q,p]
                            + F_345(r,p,q)
                            + F_345(r,q,p)
                            );
      };
      return .5*temp/r;
  },  kp-p,km  )(y); },  k, km    )(x);                 //*/
  //\
  return _1A;
  _1B =
  remap([&](double p, double pd) {                      // p=[.5km,km]
    pm = km-p; 
    pm = (fabs(pm)<1e-1*( .5*km )) ? pd : pm;           // p=km (!)
    return remap([&](double q, double qd) {             // q=[km-p,min(p,kp-p)]
      double temp = 0., r = k0-p-q;
      qm = km-q;
      temp += lga(p*q/(pm*qm))*(
                                 F_123(p,q,r)           // 1B[p,q]
                               + F_123(q,p,r)           // 1B[q,p]
                               + F_345(r,p,q)
                               + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  km-p,min(p,kp-p) )(y); }, .5*km,km   )(x);        //*/
  //\
  return _1B;
  _1D =                                                 // corner, +1C
  remap([&](double pq, double pqd) {                    // p-q=[k-min(k,km),k+min(k,km)]
    return remap([&](double r, double rd) {             // r=[|pq-k|,min(k,km)]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      pm = km-p; qm = km-q;
      pp = kp-p; qp = kp-q;
      double l=k0-p;
      double v=k0-q;
      //\
      Padding for log singularities:
      double h = fabs( min(k,km)-fabs(pq-k) );          // p=kp,q=km (!!)
      //double h = fabs( 2.*min(k,km) );// p=kp,q=km (!!)
      if (pq>k) { pp = (fabs(pp)<1e-1*h) ? rd/2. : pp; }
      if (pq<k) { qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; }
      double lg1 = lga( pp*qp/(p*q) ),
             lg2 = lga( pp*v/(qm*p) ),
             lg3 = lga( p*q/( l*v ) );

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q-MOT2,s2))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 1D[p,q]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q-MOT1,s1))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 1C[q,p]
      temp += pow(v/k0,n)*F_14(l,p)*( (.5+f(q-MOT5,s5))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); // from (*)
      temp += pow(v/k0,m)*F_25(l,p)*( (.5+f(q-MOT4,s4))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); //      (**)
      // Adjust logs:
      lg1 += lg3; lg2 += lg3;

      temp -= pow(v/k0,n)*F_14(l,p)*( (.5+f(v-MOT2,s2))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
      // (***)                      - (.5+f(q,s5))*lga( pm*qm/(p*q) )
                                    ); // 4C[l,v]
      temp -= pow(v/k0,m)*F_25(l,p)*( (.5+f(v-MOT1,s1))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
      // (****)                     - (.5+f(q,s4))*lga( pm*qm/(p*q) )
                                    ); // 4B[v,l]
      temp -= pow(q/k0,n)*F_14(p,l)*( (.5+f(v-MOT5,s5))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
                                    ); //      (***)
      temp -= pow(q/k0,m)*F_25(p,l)*( (.5+f(v-MOT4,s4))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
                                    ); //      (****)
      return .25*temp/r;
  },  fabs(pq-k),min(k,km)    )(y);
  }, k-min(k,km),k+min(k,km)  )(x);                     //*/
  //\
  return _1D;

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
  remap([&](double p, double pd) {                      // p=[kp,+inf)
    pm = km-p; pp = kp-p;
    pp = (fabs(pp)<.5*kp) ? pd : pp;                    // p=kp (!)
    return remap([&](double q, double qd) {             // q=(-inf,km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      temp += lga(pm*qm/(pp*qp))*(
                                   F_123(p,q,r)         // 2A[p,q]
                                 + F_123(q,p,r)         // 3A[q,p]
                                 + F_345(r,p,q)
                                 + F_345(r,q,p)
                                 );
      return .5*temp/r;
  },  km-p  )(y); },  kp )(x);                          //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) {                      // p=[kp,+inf)
    pp = kp-p;
    pp = (fabs(pp)<.5*kp) ? pd : pp;                    // p=kp (!)
    return remap([&](double q, double qd) {             // q=[km-p,kp-p]
      double temp = 0., r = k0-p-q;
      qp = kp-q;
      temp += lga(p*q/(pp*qp))*(
                                 F_123(p,q,r)           // 2B[p,q]
                               + F_123(q,p,r)           // 3B[q,p]
                               + F_345(r,p,q)
                               + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  km-p,kp-p  )(y); },  kp )(x);                     //*/
  //\
  return _2B;
  _2C =
  remap([&](double p, double pd) {                      // p=[km,kp]
    pm = km-p; pp = kp-p;
    pm = (fabs(pm)<1e-1*k) ? pd : pm;                   // p=km (!)
    pp = (fabs(pp)<1e-1*k) ? pd : pp;                   // p=kp (!)
    return remap([&](double q, double qd) {             // q=(-inf,km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      double lg1 = lga( pm*qm/(p*q) ),
             lg2 = lga( pm*v/(qp*p) );

       temp += pow(q/k0,n)*F_14(p,l)*( (-f(fabs(q)-sgn(q)*MOT2,s2))*lg1
                                     + (+f(fabs(r)-sgn(r)*MOT3,s3))*lg2
       // (*)                        - (.5+f(v,s5))*lga( pp*qp/(l*v) )
                                     ); // 2C[p,q]
       temp += pow(q/k0,m)*F_25(p,l)*( (-f(fabs(q)-sgn(q)*MOT1,s1))*lg1
                                     + (+f(fabs(r)-sgn(r)*MOT3,s3))*lg2
       // (**)                       - (.5+f(v,s4))*lga( pp*qp/(l*v) )
                                     ); // 3C[q,p]
       temp += pow(v/k0,n)*F_14(l,p)*( (-f(fabs(q)-sgn(q)*MOT5,s5))*lg1
                                     + (+f(fabs(r)-sgn(r)*MOT3,s3))*lg2
                                     ); // from (*)
       temp += pow(v/k0,m)*F_25(l,p)*( (-f(fabs(q)-sgn(q)*MOT4,s4))*lg1
                                     + (+f(fabs(r)-sgn(r)*MOT3,s3))*lg2
                                     ); //      (**)

      //temp -= pow(q/k0,n)*F_14(p,l)*.5*L_div(q,v,qp,qm);
      //temp -= pow(q/k0,m)*F_25(p,l)*.5*L_div(q,v,qp,qm);
      //temp -= pow(v/k0,n)*F_14(l,p)*.5*L_div(v,q,qm,qp);
      //temp -= pow(v/k0,m)*F_25(l,p)*.5*L_div(v,q,qm,qp);
      return .5*temp/r;
  },  km-p  )(y)

// This next segment is very messy :-(
//
  + remap([&](double q, double qd) {             // q=[-E*kp,km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp -= pow(q/k0,n)*F_14(p,l)*.5*L_div(q,v,qp,qm);
      temp -= pow(q/k0,m)*F_25(p,l)*.5*L_div(q,v,qp,qm);
      return .5*temp/r;
  },  -E*kp,km-p  )(y)
  + remap([&](double q, double qd) {             // q=(-inf,-E*kp]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp -= pow(q/k0,n)*F_14(p,l)*.5*L_div(q,v,qp,qm); // [A]
      temp -= pow(q/k0,m)*F_25(p,l)*.5*L_div(q,v,qp,qm); // [B]
      return .5*temp/r;
  },  -E*kp  )(y)
  + remap([&](double q, double qd) {             // q=[k0-E*kp,km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp -= pow(v/k0,n)*F_14(l,p)*.5*L_div(v,q,qm,qp);
      temp -= pow(v/k0,m)*F_25(l,p)*.5*L_div(v,q,qm,qp);
      return .5*temp/r;
  },  k0-E*kp,km-p  )(y)
  + remap([&](double q, double qd) {             // q=(-inf,k0-E*kp]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp -= pow(v/k0,n)*F_14(l,p)*.5*L_div(v,q,qm,qp); // [C]
      temp -= pow(v/k0,m)*F_25(l,p)*.5*L_div(v,q,qm,qp); // [D]
      return .5*temp/r;
  },  k0-E*kp  )(y)

    -   (
        F_14(p,k0-p)*L_int(n,k0-p,-E*kp)+ // A
        F_25(p,k0-p)*L_int(m,k0-p,-E*kp)+ // B
        F_14(k0-p,p)*L_int(n,p,E*kp)+     // C
        F_25(k0-p,p)*L_int(m,p,E*kp)      // D
        )
    ; },  km,kp )(x);                       //*/
  //\
  return _2C;
  _2D =
  remap([&](double p, double pd) {                      // p=[km,kp]
    return remap([&](double q, double qd) {             // q=[km-p,0]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      //double q_=q;
      double q_ = (fabs(q)<1e-1*fabs(km-p)) ? qd : q;   // q=0 (!)
      double l=k0-p, p2=p*p;
      double v=k0-q, q2=q*q;

      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r-MOT3,s3))*lga( q_*v/( qp*qm ) )
                                        ); // 2D[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r-MOT3,s3))*lga( q_*v/( qp*qm ) ) 
                                        ); // 3D[p,q]

      //temp += pow(q/k0,m)*F_25(p,l)*( (m==2) ? K2/(8.*q*q) : 0.); // Case: m=2
      //temp += pow(v/k0,m)*F_25(l,p)*( (m==2) ? K2/(8.*v*v) : 0.);
      return .5*temp/r;
  },  km-p, 0.  )(y); },  km,kp    )(x);                //*/
  //\
  return _2D;

  //  #] from '2','3'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) {                      // p=[kp,+inf)
    pm = km-p; pp = kp-p;
    pp = (fabs(pp)<.5*kp) ? pd : pp;                    // p=kp (!)
    return remap([&](double q, double qd) {             // q=[p,+inf)
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      temp += lga(pm/(pp))*(
                                   F_123(p,q,r)         // 4A[p,q]
                                 + F_123(q,p,r)         // 4A[q,p]
                                 + F_345(r,p,q)
                                 + F_345(r,q,p)
                                 );
      return .5*temp/r;
  },  kp  )(y); },  kp    )(x);                         //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) {                    // q=[km,kp]
    qm = km-q; qp = kp-q;
    qm = (fabs(qm)<1e-1*k) ? qd : qm;                 // q=km (!)
    qp = (fabs(qp)<1e-1*k) ? qd : qp;                 //  =kp (!)
    return remap([&](double p, double pd) {           // p=[k0,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l = k0-p, v=k0-q;
      //l = (fabs(l)<1e-1*k0) ? pd : l;
      double lg1 = lga( pm*qm/(p*q) ),
             lg2 = lga( qm*l/(pp*q) );

      temp += pow(p/k0,n)*F_14(q,v)*( (+f(fabs(p)-sgn(p)*MOT2,s2))*lg1
                                    + (-f(fabs(r)-sgn(r)*MOT3,s3))*lg2
      // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                    ); // 4C[q,p]
      temp += pow(p/k0,m)*F_25(q,v)*( (+f(fabs(p)-sgn(p)*MOT1,s1))*lg1
                                    + (-f(fabs(r)-sgn(r)*MOT3,s3))*lg2
      // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                    ); // 4B[p,q]
      temp += pow(l/k0,n)*F_14(v,q)*( (+f(fabs(p)-sgn(p)*MOT5,s5))*lg1
                                    + (-f(fabs(r)-sgn(r)*MOT3,s3))*lg2
                                    ); // from (*)
      temp += pow(l/k0,m)*F_25(v,q)*( (+f(fabs(p)-sgn(p)*MOT4,s4))*lg1
                                    + (-f(fabs(r)-sgn(r)*MOT3,s3))*lg2
                                    ); //      (**)

      //temp += pow(p/k0,n)*F_14(q,v)*.5*L_div(p,l,pp,pm);
      //temp += pow(p/k0,m)*F_25(q,v)*.5*L_div(p,l,pp,pm);
      //temp += pow(l/k0,n)*F_14(v,q)*.5*L_div(l,p,pm,pp);
      //temp += pow(l/k0,m)*F_25(v,q)*.5*L_div(l,p,pm,pp);

      return .5*temp/r;
  },  k0  )(y)

    + remap([&](double p, double pd) {           // p=[k0,E*kp]
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l = k0-p, v=k0-q;
      temp += pow(p/k0,n)*F_14(q,v)*.5*L_div(p,l,pp,pm);
      temp += pow(p/k0,m)*F_25(q,v)*.5*L_div(p,l,pp,pm);
      return .5*temp/r;
    },  k0, E*kp  )(y)
    + remap([&](double p, double pd) {           // p=[E*kp,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l = k0-p, v=k0-q;
      temp += pow(p/k0,n)*F_14(q,v)*.5*L_div(p,l,pp,pm);
      temp += pow(p/k0,m)*F_25(q,v)*.5*L_div(p,l,pp,pm);
      return .5*temp/r;
    },  E*kp  )(y)
    + remap([&](double p, double pd) {           // p=[k0,E*kp+k0]
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l = k0-p, v=k0-q;
      temp += pow(l/k0,n)*F_14(v,q)*.5*L_div(l,p,pm,pp);
      temp += pow(l/k0,m)*F_25(v,q)*.5*L_div(l,p,pm,pp);
      return .5*temp/r;
    },  k0, E*kp+k0  )(y)
    + remap([&](double p, double pd) {           // p=[E*kp+k0,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l = k0-p, v=k0-q;
      temp += pow(l/k0,n)*F_14(v,q)*.5*L_div(l,p,pm,pp);
      temp += pow(l/k0,m)*F_25(v,q)*.5*L_div(l,p,pm,pp);
      return .5*temp/r;
    },  E*kp+k0  )(y)

    - (
        F_14(q,k0-q)*L_int(n,k0-q,E*kp)+
        F_25(q,k0-q)*L_int(m,k0-q,E*kp)+
        F_14(k0-q,q)*L_int(n,q,-E*kp)+
        F_25(k0-q,q)*L_int(m,q,-E*kp)
        )
    ; }
  ,  km,kp    )(x);                                   //*/
  //\
  return _4B;
  if (k0>3.*k) {
    //
    _4B+=
    remap([&](double p, double pd) {                  // q=[km,kp]
      return remap([&](double q, double qd) {         // p=[k0+km-q,k0]
        double temp = 0., r = k0-p-q;
        pm = km-p; pp = kp-p;
        qm = km-q;
        //qm = (fabs(qm)<1e-1*k) ? qd : qm;           // q=km (!)
        double l=k0-p; //l = (fabs(l)<1e-1*fabs(k)) ? pd : l;
        double v=k0-q;
        double lg1 = lga( pm*qm/(p*q) ),
               lg2 = lga( qm*l/(pp*q) );

        temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
        // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                      ); // 4C[q,p]
        temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
        // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                      ); // 4B[p,q]
        temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p-MOT5,s5))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
                                      ); // from (*)
        temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p-MOT4,s4))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
                                      ); //      (**)

        temp += ( pow(l/k0,n)*F_14(v,q)
                + pow(p/k0,n)*F_14(q,v) )*( (.5+f(r-MOT3,s3))*(lg2-lg1) // log(l*p/(pp*pm))
                                          ); // 2D[v,l]
        temp += ( pow(l/k0,m)*F_25(v,q)
                + pow(p/k0,m)*F_25(q,v) )*( (.5+f(r-MOT3,s3))*(lg2-lg1)
                                          ); // 3D[l,v]
        return .5*temp/r;
    },  k0+km-p,kp  )(y); },  2.*km,k0    )(x); //*/
    //
    _4B+=
    remap([&](double p, double pd) {                  // q=[km,kp]
      return remap([&](double q, double qd) {         // p=[2kp-q,k0+km-q]
        double temp = 0., r = k0-p-q;
        pm = km-p; pp = kp-p;
        qm = km-q;
        double h = fabs( min(kp,k0+km-p) - max(km,2.*kp-p) );
        if (p>2.*kp-km) { qm = (fabs(qm)<1e-1*h) ? qd : qm; }// q=km (!)
        double l = k0-p, v=k0-q;
        l  = ( fabs(l)<1e-1*km) ? pd : l;
        pp = (fabs(pp)<1e-1*km) ? pd : pp;
        double lg1 = lga( pm*qm/(p*q) ),
               lg2 = lga( qm*l/(pp*q) ),
               lg3 = lga( p*q/( l*v ) );

        temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
        // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                      ); // 4C[q,p]
        temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
        // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                      ); // 4B[p,q]
        temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p-MOT5,s5))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
                                      ); // from (*)
        temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p-MOT4,s4))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
                                      ); //      (**)
        lg1 += lg3; lg2 += lg3;

        temp -= pow(l/k0,n)*F_14(v,q)*( (.5+f(l-MOT2,s2))*lg1
                                      - (.5+f(r-MOT3,s3))*lg2
        // (***)                      - (.5+f(p,s5))*lga( pp*qp/(q*p) )
                                      ); // 1D[v,l]
        temp -= pow(l/k0,m)*F_25(v,q)*( (.5+f(l-MOT1,s1))*lg1
                                      - (.5+f(r-MOT3,s3))*lg2
        // (****)                     - (.5+f(p,s4))*lga( pp*qp/(q*p) )
                                      ); // 1C[l,v]
        temp -= pow(p/k0,n)*F_14(q,v)*( (.5+f(l-MOT5,s5))*lg1
                                      - (.5+f(r-MOT3,s3))*lg2
                                      ); //      (***)
        temp -= pow(p/k0,m)*F_25(q,v)*( (.5+f(l-MOT4,s4))*lg1
                                      - (.5+f(r-MOT3,s3))*lg2
                                      ); //      (****)
        return .5*temp/r;
    //},  2.*kp-q,k0+km-q  )(y); },  km,kp    )(x); //*/
    },  max(km,2.*kp-p),min(kp,k0+km-p) )(y); },  kp,k0    )(x); //*/
    //
  } else if (k0<3.*k) {
    //
    _4B+=
    remap([&](double p, double pd) {                  // p=[kp,k0]
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<1e-1*fabs(k0-kp)) ? pd : pp;     // p=kp (!)
      double l = k0-p; 
      l = (fabs(l)<1e-1*fabs(k0-kp)) ? pd : l;        // p=k0 (!)
      return remap([&](double q, double qd) {         // q=[k0+km-p,kp]
        double temp = 0., r = k0-p-q;
        qm = km-q, qp=kp-q;
        double v=k0-q;
        double lg1 = lga( pm*qm/(p*q) ),
               lg2 = lga( qm*l/(pp*q) );

        temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
        // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                      ); // 4C[q,p]
        temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
        // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                      ); // 4B[p,q]
        temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p-MOT5,s5))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
                                      ); // from (*)
        temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p-MOT4,s4))*lg1
                                      + (.5+f(r-MOT3,s3))*lg2
                                      ); //      (**)

        temp += ( pow(l/k0,n)*F_14(v,q)
                + pow(p/k0,n)*F_14(q,v) )*( (.5+f(r-MOT3,s3))*(lg2-lg1) // log(l*p/(pp*pm))
                                          ); // 2D[v,l]
        temp += ( pow(l/k0,m)*F_25(v,q)
                + pow(p/k0,m)*F_25(q,v) )*( (.5+f(r-MOT3,s3))*(lg2-lg1)
                                          ); // 3D[l,v]
        return .5*temp/r;
    },  k0+km-p, kp  )(y); },  kp,k0    )(x);         //*/
  }
  //\
  return _4B;
  _4D =
  remap([&](double pq, double pqd) {                  // p-q=[max(0,k-km),k]
    return remap([&](double r, double rd) {           // r=[p-q-k,k-p+q]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      pm = km-p; qm = km-q;
      pp = kp-p; qp = kp-q;
      double l=k0-p, p2=p*p, l2=l*l;
      double v=k0-q, q2=q*q, v2=v*v;

      double h = fabs( k-(pq) );
      qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q-MOT2,s2))*lga( qm*qp/q2 )
                                    - (.5+f(v-MOT5,s5))*lga( qm*qp/v2 ) ); // 4D[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lga( pm*pp/p2 )
                                    - (.5+f(l-MOT4,s4))*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lga( pm*pp/p2 )
                                    - (.5+f(l-MOT5,s5))*lga( pm*pp/l2 ) ); // 4D[q,p]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q-MOT1,s1))*lga( qm*qp/q2 )
                                    - (.5+f(v-MOT4,s4))*lga( qm*qp/v2 ) );

      temp -= pow(v/k0,n)*F_14(l,p)*( (.5+f(v-MOT2,s2))*lga( qm*qp/v2 )
                                    - (.5+f(q-MOT5,s5))*lga( qm*qp/q2 ) ); // 4D[l,v]
      temp -= pow(l/k0,m)*F_25(v,q)*( (.5+f(l-MOT1,s1))*lga( pm*pp/l2 )
                                    - (.5+f(p-MOT4,s4))*lga( pm*pp/p2 ) );

      temp -= pow(l/k0,n)*F_14(v,q)*( (.5+f(l-MOT2,s2))*lga( pm*pp/l2 )
                                    - (.5+f(p-MOT5,s5))*lga( pm*pp/p2 ) ); // 4D[v,l]
      temp -= pow(v/k0,m)*F_25(l,p)*( (.5+f(v-MOT1,s1))*lga( qm*qp/v2 )
                                    - (.5+f(q-MOT4,s4))*lga( qm*qp/q2 ) );
      return .25*temp/r;
  },  0.,k-(pq) )(y); }, 0.,k    )(x);                //*/
  //
  if (k0<3.*k) { 
    //
    _4D +=
    remap([&](double pq, double pqd) {                // p-q=[0,k-km]
      return remap([&](double r, double rd) {         // r=[km,k-pq]
        double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
        pm = km-p; qm = km-q;
        pp = kp-p; qp = kp-q;
        double h = fabs(k-pq-km);
        double l=k0-p, p2=p*p, l2=l*l;
        double v=k0-q, q2=q*q, v2=v*v;
        qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;          // q=km (!)

        temp -= lga(qp*pp/(p*q))*(
                                   F_123(p,q,r) 
                                 + F_123(q,p,r) 
                                 + F_345(r,p,q)
                                 + F_345(r,q,p)
                                 );
        return .25*temp/r;
    },  km,k-pq  )(y); }, 0.,k-km    )(x);            //*/
    //
  }
  //\
  return _4D;

  res += 
  checker(_1A) + 
  checker(_1B) + 
  checker(_1D) + 
  checker(_2A) + 
  checker(_2B) + 
  checker(_2C) + 
  checker(_2D) + 
  checker(_4A) + 
  checker(_4B) + 
  checker(_4D);

  } else {

  double _2A, _2B,      // 3 by refl.
         _4A, _4B,      _4D,
         _5A, _5B,
         _6A, _6B, _6C; // 7 by refl.

  //  #[ from '2'

  _2A =
  remap([&](double p, double pd) {                    // p=[km,0]
    return remap([&](double q, double qd) {           // q=(-inf,km]
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double lg1 = lga( p*l/(pp*pm) );
      temp += ( pow(p/k0,n)*F_14(q,v)
              + pow(l/k0,n)*F_14(v,q) )*( (.5+f(r-MOT3,s3))*lg1
                                       ); // 2A[p,q]
      temp += ( pow(p/k0,m)*F_25(q,v)
              + pow(l/k0,m)*F_25(v,q) )*( (.5+f(r-MOT3,s3))*lg1
                                        ); // 3A[q,p]
      return .5*temp/r;
  },  km  )(y); },  km,0.    )(x);
  //
  _2A +=
  remap([&](double p, double pd) {                    // p=[0,kp]
    return remap([&](double q, double qd) {           // q=(-inf,km-p]
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double lg1 = lga( p*l/(pp*pm) );
      temp += ( pow(p/k0,n)*F_14(q,v)
              + pow(l/k0,n)*F_14(v,q) )*( (.5+f(r-MOT3,s3))*lg1
                                        ); // 2A[p,q]
      temp += ( pow(p/k0,m)*F_25(q,v)
              + pow(l/k0,m)*F_25(v,q) )*( (.5+f(r-MOT3,s3))*lg1
                                        ); // 3A[p,q]
      return .5*temp/r;
  },  km-p  )(y); },  0.,kp    )(x);                  //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) {                    // p=[0,k0]
    return remap([&](double q, double qd) {           // q=[km-p,km]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double h = fabs( p );
      qm = (fabs(qm)<1e-1*h) ? qd : qm;
      double lg1 = lga( pm*qm/(p*q) ),
             lg2 = lga( qm*l/(pp*q) );

      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 2B[p,q]
      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 3B[q,p]
      temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p-MOT4,s4))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); // from (*)
      temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p-MOT5,s5))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); //      (**)
      return .5*temp/r;
  },  km-p,km  )(y); },  0.,k0    )(x);
  //
  _2B +=
  remap([&](double p, double pd) {                    // p=[k0,kp]
    return remap([&](double q, double qd) {           // q=[km-p,k0+km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double h = fabs( kp-k0 );
      pp = (fabs(pp)<1e-1*h) ? pd : pp;
      double lg1 = lga( pm*qm/(p*q) ),
             lg2 = lga( qm*l/(pp*q) );

      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 2B[p,q]
      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 3B[q,p]
      temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p-MOT4,s4))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); // from (*)
      temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p-MOT5,s5))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); //      (**)
      return .5*temp/r;
  },  km-p,k0+km-p  )(y); },  k0,kp    )(x);
  //
  _2B +=
  remap([&](double pq, double pqd) {                  // p-q=[k+km,k-km]
    return remap([&](double r, double rd) {           // r=[|p-q-k|,-km]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      double h = fabs( -km-fabs(pq-k) );
      if (pq>k) { pp = (fabs(pp)<1e-1*h) ? rd/2. : pp; }
      if (pq<k) { qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; }
      double l=k0-p, v=k0-q;
      double lg1 = lga( pm*qm/(p*q) ),
             lg2 = lga( qm*l/(pp*q) ),
             lg3 = lga( p*q/( l*v ) );

      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p-MOT1,s1))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 2B[p,q]
      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p-MOT2,s2))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 3B[q,p]
      temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p-MOT4,s4))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); // from (*)
      temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p-MOT5,s5))*lg1
                                    + (.5+f(r-MOT3,s3))*lg2
                                    ); //      (**)

      lg1 += lg3; lg2 += lg3;

      temp -= pow(l/k0,n)*F_14(v,q)*( (.5+f(l-MOT2,s2))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
                                    ); // 4B'[v,l]
      temp -= pow(l/k0,m)*F_25(v,q)*( (.5+f(l-MOT1,s1))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
                                    ); // 4C'[l,v]
      temp -= pow(p/k0,n)*F_14(q,v)*( (.5+f(l-MOT5,s5))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
                                    );
      temp -= pow(p/k0,m)*F_25(q,v)*( (.5+f(l-MOT4,s4))*lg1
                                    - (.5+f(r-MOT3,s3))*lg2
                                    );
      return .25*temp/r;
  },  fabs(pq-k),-km  )(y); }, k+km,k-km  )(x);       //*/
  //\
  return _2B;

  //  #] from '2'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) {                    // p=[kp,+inf)
    return remap([&](double q, double qd) {           // q=[kp,p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      pp = kp-p; pm = km-p;
      /* not so confident about this way:
      double l=k0-p, v=k0-q;
      temp += pow(q/k0,n)*F_14(p,l)*( f(q,s2)-f(r,s3) )*lga(qp/q) ;
      temp += pow(v/k0,n)*F_14(l,p)*( f(q,s5)-f(r,s3) )*lga(qp/q) ;
      temp += pow(p/k0,m)*F_25(q,v)*( f(p,s1)-f(r,s3) )*lga(pp/p) ;
      temp += pow(l/k0,m)*F_25(v,q)*( f(p,s4)-f(r,s3) )*lga(pp/p) ;*/

      double p2=p*p, q2=q*q,
             l=k0-p, l2=l*l,
             v=k0-q, v2=v*v;

      temp += pow(q/k0,n)*F_14(p,l)*( f(fabs(q)-sgn(q)*MOT2,s2)*lga( qm*qp/q2 )
                                    + f(fabs(v)-sgn(v)*MOT5,s5)*lga( qm*qp/v2 ) ); // 4A[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( f(fabs(p)-sgn(p)*MOT1,s1)*lga( pm*pp/p2 )
                                    + f(fabs(l)-sgn(l)*MOT4,s4)*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( f(fabs(p)-sgn(p)*MOT2,s2)*lga( pm*pp/p2 )
                                    + f(fabs(l)-sgn(l)*MOT5,s5)*lga( pm*pp/l2 ) ); // 4A[q,p]
      temp += pow(q/k0,m)*F_25(p,l)*( f(fabs(q)-sgn(q)*MOT1,s1)*lga( qm*qp/q2 )
                                    + f(fabs(v)-sgn(v)*MOT4,s4)*lga( qm*qp/v2 ) );

      //temp += pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm); // q>E*kp to be re-included
      //temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm); // p>E*kp "
      //temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
      //temp += pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);

      temp += lga(qm*pm/( p*q ))*(
                                  F_123(p,q,r)
                                + F_123(q,p,r)
                                + F_345(r,p,q)
                                + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  kp,p  )(y); },  kp    )(x);                     //*/

// This next segment is very messy :-(
//
   _4A += remap([&](double q, double qd) {             // q=[kp,+inf)
    return 
      remap([&](double p, double pd) {          // p=[kp,E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        return .5*temp/r;
      }, kp, E*kp )(y)+
      remap([&](double p, double pd) {                // p=[E*kp,+inf)
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        return .5*temp/r; 
      }, E*kp )(y)
      - F_25(q,k0-q)*L_int(m,k0-q,E*kp)*2.
      - F_14(q,k0-q)*L_int(n,k0-q,E*kp)*2.;
  }, kp    )(x);                     //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) {                    // q=[k0+km,kp]
    return remap([&](double p, double pd) {           // p=[kp,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r-MOT3,s3))*lga( q*v/( qp*qm ) )
                                        ); // 4B[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r-MOT3,s3))*lga( q*v/( qp*qm ) ) 
                                        ); // 4C[q,p]
      return .5*temp/r;
  },  kp  )(y); },  k0+km,kp    )(x);
  //
  _4B +=
  remap([&](double q, double qd) {                    // q=[0,k0+km]
    return remap([&](double p, double pd) {           // p=[kp,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r-MOT3,s3))*lga( q*v/( qp*qm ) )
                                        ); // 4B[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r-MOT3,s3))*lga( q*v/( qp*qm ) ) 
                                        ); // 4C[q,p]
      return .5*temp/r;
  },  kp  )(y); },  0.,k0+km    )(x);
  //
  _4B +=
  remap([&](double q, double qd) {                    // q=[km,0]
    return remap([&](double p, double pd) {           // p=[kp-q,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r-MOT3,s3))*lga( q*v/( qp*qm ) )
                                        ); // 4B[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r-MOT3,s3))*lga( q*v/( qp*qm ) ) 
                                        ); // 4C[q,p]
      return .5*temp/r;
  },  kp-q  )(y); }, km, 0.  )(x);                    //*/
  //\
  return _4B;
  //
  _4D =
  remap([&](double q, double qd) {                    // q=[0,kp]
    return remap([&](double p, double pd) {           // p=[kp-q,kp]
      double temp = 0., r = k0-p-q;
      qp = kp-q;
      pp = kp-p;
      //\
      Padding for log singularities:
      double h = fabs(q );
      // if (p>q) 
      pp = (fabs(pp)<1e-1*h) ? pd : pp;
      //if (p<q) { qp = (fabs(qp)<1e-1*h) ? qd : qm; }

      temp += lga( pp*qp/(p*q) )*( F_123(p,q,r) + F_345(r,p,q) ); // 4D[p,q]
      return .5*temp/r;
  },  kp-q,kp  )(y); },  0.,kp    )(x);               //*/
  //\
  return _4D;

  //  #] from '4'\
      #[ from '5'

  _5A =
  remap([&](double q, double qd) {                    // q=(-inf,km]
    return remap([&](double p, double pd) {           // p=[q,km]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      /* not so confident about this way:
      double l=k0-p, v=k0-q;
      temp += pow(q/k0,n)*F_14(p,l)*( f(q,s2)-f(r,s3) )*lga(qm/q) ;
      temp += pow(v/k0,n)*F_14(l,p)*( f(q,s5)-f(r,s3) )*lga(qm/q) ;
      temp += pow(p/k0,m)*F_25(q,v)*( f(p,s1)-f(r,s3) )*lga(pm/p) ;
      temp += pow(l/k0,m)*F_25(v,q)*( f(p,s4)-f(r,s3) )*lga(pm/p) ; */
      double p2=p*p, q2=q*q,
             l=k0-p, l2=l*l,
             v=k0-q, v2=v*v;

      temp -= pow(q/k0,n)*F_14(p,l)*( f(fabs(q)-sgn(q)*MOT2,s2)*lga( qm*qp/q2 )
                                    + f(fabs(v)-sgn(v)*MOT5,s5)*lga( qm*qp/v2 ) ); // 6B[p,q]
      temp -= pow(p/k0,m)*F_25(q,v)*( f(fabs(p)-sgn(p)*MOT1,s1)*lga( pm*pp/p2 )
                                    + f(fabs(l)-sgn(l)*MOT4,s4)*lga( pm*pp/l2 ) );

      temp -= pow(p/k0,n)*F_14(q,v)*( f(fabs(p)-sgn(p)*MOT2,s2)*lga( pm*pp/p2 )
                                    + f(fabs(l)-sgn(l)*MOT5,s5)*lga( pm*pp/l2 ) ); // 7B[q,p]
      temp -= pow(q/k0,m)*F_25(p,l)*( f(fabs(q)-sgn(q)*MOT1,s1)*lga( qm*qp/q2 )
                                    + f(fabs(v)-sgn(v)*MOT4,s4)*lga( qm*qp/v2 ) );

      //temp -= pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm);
      //temp -= pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
      //temp -= pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
      //temp -= pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);

      temp += lga(qp*pp/(p*q))*(
                                  F_123(p,q,r)
                                + F_123(q,p,r)
                                + F_345(r,p,q)
                                + F_345(r,q,p)
                               );

      return .5*temp/r;
  },  q,km  )(y); },  km    )(x);                     //*/
  _5A += remap([&](double q, double qd) {             // q=(-inf,km]
    return 
      remap([&](double p, double pd) {                // p=[-E*kp,km]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        temp -= pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp -= pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        return .5*temp/r;
      }, -E*kp, km )(y)+
      remap([&](double p, double pd) {                // p=(-inf,-E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        temp -= pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp -= pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        return .5*temp/r; 
      }, -E*kp )(y)
      - F_25(q,k0-q)*L_int(m,k0-q,-E*kp)*2.
      - F_14(q,k0-q)*L_int(n,k0-q,-E*kp)*2.;
  }, km    )(x);                     //*/
  //\
  return _5A;
  //
  _5B =
  remap([&](double q, double qd) {                    // q=[km,0]
    return remap([&](double p, double pd) {           // p=[km,km-q]
      double temp = 0., r = k0-p-q;
      pm = km-p;
      qm = km-q;
      double h = fabs(q );
      pm = (fabs(pm)<1e-1*h) ? pd : pm;               // p=km (!)

      temp += lga( pm*qm/(p*q) )*( F_123(p,q,r) + F_345(r,p,q) );   // 5B[p,q]

      return .5*temp/r;
  },  km,km-q  )(y); },  km,0.    )(x);               //*/

  //  #] from '5'\
      #[ from '6'

  _6A =
  remap([&](double p, double pd) {                    // p=[k,+inf)
    return remap([&](double q, double qd) {           // q=[kp-p,km]
      double temp = 0., r = k0-p-q;
      pp = kp-p;
      qp = kp-q;
      double p2=p*p, q2=q*q,
             l=k0-p, l2=l*l,
             v=k0-q, v2=v*v;

      temp += lga(qp*pp/(p*q))*(
                                 F_123(p,q,r)
                               + F_123(q,p,r)
                               + F_345(r,p,q)
                               + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  kp-p,km  )(y); },  k    )(x);                   //*/
  //\
  return _6A;
  //
  _6B =
  remap([&](double pq, double pqd) {                  // p-q=[k,+inf)
    return remap([&](double r, double rd) {           // r=[k-pq,pq-k]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      double p2=p*p, q2=q*q,
             l=k0-p, l2=l*l,
             v=k0-q, v2=v*v;

      double h = 2.*fabs( pq-k );
      if (r>0) { pp = (fabs(pp)<1e-1*h) ? rd/2. : pp; }
      if (r<0) { qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; }

      temp -= pow(q/k0,n)*F_14(p,l)*( f(fabs(q)-sgn(q)*MOT2,s2)*lga( qm*qp/q2 )
                                    + f(fabs(v)-sgn(v)*MOT5,s5)*lga( qm*qp/v2 ) ); // 6B[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( f(fabs(p)-sgn(p)*MOT1,s1)*lga( pm*pp/p2 )
                                    + f(fabs(l)-sgn(l)*MOT4,s4)*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( f(fabs(p)-sgn(p)*MOT2,s2)*lga( pm*pp/p2 )
                                    + f(fabs(l)-sgn(l)*MOT5,s5)*lga( pm*pp/l2 ) ); // 7B[q,p]
      temp -= pow(q/k0,m)*F_25(p,l)*( f(fabs(q)-sgn(q)*MOT1,s1)*lga( qm*qp/q2 )
                                    + f(fabs(v)-sgn(v)*MOT4,s4)*lga( qm*qp/v2 ) );

      //temp -= pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm);
      //temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
      //temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
      //temp -= pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);

      return .25*temp/r;
  },  k-pq,pq-k  )(y); }, k    )(x);                  //*/
  _6B += remap([&](double q, double qd) {             // q=(-inf,km]
     return 
      remap([&](double p, double pd) {                // p=[kp,E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,m)*F_25(l,p)*L_div(v,q,qm,qp);
        temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,n)*F_14(l,p)*L_div(v,q,qm,qp);
        return .125*temp/r;
      }, kp,E*kp )(y)+
      remap([&](double p, double pd) {                // p=[E*kp,+inf)
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,m)*F_25(l,p)*L_div(v,q,qm,qp);
        temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,n)*F_14(l,p)*L_div(v,q,qm,qp);
        return .125*temp/r; 
      }, E*kp )(y)
      - F_25(q,k0-q)*L_int(m,k0-q,E*kp)
      - F_14(q,k0-q)*L_int(n,k0-q,E*kp);
  }, -E*kp,km    )(x);                     //
  _6B += remap([&](double q, double qd) {             // q=(-inf,km]
     return 
      remap([&](double p, double pd) {                // p=[kp,E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,m)*F_25(l,p)*L_div(v,q,qm,qp);
        temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,n)*F_14(l,p)*L_div(v,q,qm,qp);
        return .125*temp/r;
      }, kp,E*kp )(y)+
      remap([&](double p, double pd) {                // p=[kp,E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp += pow(p/k0,m)*F_25(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,m)*F_25(l,p)*L_div(v,q,qm,qp);
        temp += pow(p/k0,n)*F_14(q,v)*L_div(p,l,pp,pm);
        temp -= pow(v/k0,n)*F_14(l,p)*L_div(v,q,qm,qp);
        return .125*temp/r;
      }, E*kp )(y)
      - F_25(q,k0-q)*L_int(m,k0-q,E*kp)
      - F_14(q,k0-q)*L_int(n,k0-q,E*kp);
  }, -E*kp    )(x);                     //*/
  _6B += remap([&](double p, double pd) {             // p=[kp,E*kp]
     return 
      remap([&](double q, double qd) {                // q=[-E*kp,km]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp -= pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,n)*F_14(v,q)*L_div(l,p,pm,pp);
        temp -= pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,m)*F_25(v,q)*L_div(l,p,pm,pp);
        return .125*temp/r;
      }, -E*kp,km )(y)+
      remap([&](double q, double qd) {                // q=(-inf,-E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp -= pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,n)*F_14(v,q)*L_div(l,p,pm,pp);
        temp -= pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,m)*F_25(v,q)*L_div(l,p,pm,pp);
        return .125*temp/r; 
      }, -E*kp )(y)
      - F_25(p,k0-p)*L_int(m,k0-p,-E*kp)
      - F_14(p,k0-p)*L_int(n,k0-p,-E*kp);
  }, kp,E*kp    )(x);                     //*/
  _6B += remap([&](double p, double pd) {             // p=[E*kp,+inf)
     return 
      remap([&](double q, double qd) {                // q=[-E*kp,km]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp -= pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,n)*F_14(v,q)*L_div(l,p,pm,pp);
        temp -= pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,m)*F_25(v,q)*L_div(l,p,pm,pp);
        return .25*temp/r;
      }, -E*kp,km )(y)+
      remap([&](double q, double qd) {                // q=(-inf,E*kp]
        double temp = 0., r=k0-p-q, v=k0-q, l=k0-p;
        pm = km-p; pp = kp-p;
        qm = km-q; qp = kp-q;
        temp -= pow(q/k0,n)*F_14(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,n)*F_14(v,q)*L_div(l,p,pm,pp);
        temp -= pow(q/k0,m)*F_25(p,l)*L_div(q,v,qp,qm);
        temp += pow(l/k0,m)*F_25(v,q)*L_div(l,p,pm,pp);
        return .25*temp/r;
      }, -E*kp )(y)
      - F_25(p,k0-p)*L_int(m,k0-p,-E*kp)*2. // NB, see also 5A
      - F_14(p,k0-p)*L_int(n,k0-p,-E*kp)*2.;
  }, E*kp    )(x);                     //*/
  /*_6B += remap([&](double q, double qd) {             // q=[kp,+inf)
    double temp = 0., v=k0-q;
    temp += (
             -F_25(q,v)*L_int(m,v,+E*kp)
             -F_14(q,v)*L_int(n,v,+E*kp)
            );
    return temp;
  }, km    )(x);                     //*/
  /*_6B += remap([&](double p, double pd) {             // p=[kp,+inf)
    double temp = 0., l=k0-p;
    temp += (
             -F_14(p,l)*L_int(n,l,-E*kp)
             -F_25(p,l)*L_int(m,l,-E*kp)
            );
    return temp;
  }, kp    )(x);                     //*/
  //\
  return _6B;
  _6C =
  remap([&](double q, double qd) {                    // q=(-inf,-k]
    return remap([&](double p, double pd) {           // p=[kp,km-q]
      double temp = 0., r = k0-p-q;
      qm = km-q;
      pm = km-p;
      double l=k0-p, v=k0-q;

      temp += lga(qm*pm/(p*q))*(
                                  F_123(p,q,r)
                                + F_123(q,p,r)
                                + F_345(r,p,q)
                                + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  kp,km-q  )(y); },  -k    )(x);                  //*/

  //  #] from '6'

  res +=
  checker(_2A) + 
  checker(_2B) + 
  checker(_4A) + 
  checker(_4B) + 
  checker(_4D) + 
  checker(_5A) + 
  checker(_5B) + 
  checker(_6A) + 
  checker(_6B) + 
  checker(_6C);
  }

  // still don't have a good way to cater for NaNs
  //if ( isinf(res)||isnan(res) ) { return -1e6;}//for inspection
  //if ( isinf(res)||isnan(res) ) { return 0.;}
  //else 
  if ( m==2 && n==0 ) {
    return .5*SQR(k0)*(km/kp)*res/k;
  }
  if ( (m==1)&&(n==1) ) {
    double temp = 0.;
    temp += ( (k0>k) ? 2.*k : 0. );
    double ep = exp(-fabs(kp)), em = exp(-fabs(km));
    temp += lga( // thermal coefficient of `vacuum'-part
         (1.-s1*ep)*(1.-s4*ep)*(1.-s2*ep)*(1.-s5*ep)/(
         (1.-s1*em)*(1.-s4*em)*(1.-s2*em)*(1.-s5*em) ));
    return .5*SQR(k0)*(km/kp)*( res - (this->OPE).T0*temp*K2/SQR(k0) )/k;
  }

  { return .5*K2*(km/kp)*res/k; }
  //else { return res/k; }
}
