#include "core.hh"
#include "quad.hh"
#include "map.hh"

using namespace std;

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
}//*/

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
// evaluation step, including the OPE

double rho11111::eval() 
{
  double res, err;

  double a1=I(0,(this->s)[1]), b1=I(2,(this->s)[1]), // tadpole ints
         a2=I(0,(this->s)[2]), b2=I(2,(this->s)[2]),
         a3=I(0,(this->s)[3]), b3=I(2,(this->s)[3]),
         a4=I(0,(this->s)[4]), b4=I(2,(this->s)[4]),
         a5=I(0,(this->s)[5]), b5=I(2,(this->s)[5]);

  if ( m==0 && n==0 ) { // (0)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 = -( a1+a2+2.*a3+a4+a5 )*.25*OOFP;
    (this->OPE).T4 = -( (b1+b2+b4+b5)*11.+b3*6.)*(k0*k0+k*k/3.)/SQR(K2)/6.*OOFP;
  } else
  if ( m==1 && n==1 ) { // (1,1)
    (this->OPE).T0 = K2/16.;
    (this->OPE).T2 = 0.;
    (this->OPE).T4 = -(( (b1+b2+b4+b5)*3.-b3*2.      )*.5*K2-
                       ( (b1+b2)*9.+b3*2.+(b4+b5)*5. )*k0*k0 )/SQR(K2)/12.*OOFP;
  } else
  if ( m==2 && n==0 ) { // (2)
    (this->OPE).T0 = -.5*K2*( +11./16. -.5 );
    (this->OPE).T2 = -( (a2+a3+a4)*k0*k0-(a2+a5)*K2*.25 )*.25/K2*OOFP;
    (this->OPE).T4 = -(((b1+b4)*3.+(b2+b3+b5)*2.)*.5*K2-
                       (b2*7.+b3+b4*9.+b5*2.)*k0*k0+
                       ((b2+b4)*11.+b3 )*k0*k0*(k0*k0+k*k/3.)/K2  )/SQR(K2)/6.*OOFP;
  } else {
    cerr << "Case: (m,n)=("<< m << ","<<n<<") out of bounds!\n";
    return 0.;
  }
  //outer f1;
  //f1.f2.R = this;
  //integrate<outer> I(f1); // do the x-integral
  //res = go(I);
  double epsabs = 1e-2, epsrel = 1e-2;
  size_t limit = 1e6;

  quad wsp1(limit);
  quad wsp2(limit);

  auto outer = make_gsl_function( [&](double x) 
  {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
          return (this->integrand)(x,y);
        } );
    gsl_integration_qag( inner, .0+1e-10,1., epsabs, 1e-3,
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag( outer, .0+1e-10,1., epsabs, 1e-2,
                       limit, 6, wsp2, &res, &err  );//*/

  double temp = 0.;
  if ( (m==2)&&(n==0) ) {
    temp += ( (k0>k) ? 1. : 0. );
    double ep = exp(-fabs(kp)), em = exp(-fabs(km));
    temp += lga( 
         (1.-((double)(this->s)[1])*ep)*(1.-((double)(this->s)[2])*ep)/
        ((1.-((double)(this->s)[1])*em)*(1.-((double)(this->s)[2])*em)) )/k;
    temp *= 1.;
  };
  return ( (  res*pow(k0,m+n)*pow(K2,-(m+n)/2.) ) )*CUBE(OOFP);
  /*return ( (  res*pow(k0,m+n)*pow(K2,-(m+n)/2.)
         //+ k0*k0*temp
        -(this->OPE).T0
      )*CUBE(OOFP) - (this->OPE).T2 )/(this->OPE).T4;*/
}

/*--------------------------------------------------------------------*/
// integrand manipulations

inline double L_div(double p, double l, double pp, double pm) {
  //double pp = kp-p, pm = km-p, l=k0-p;
  double res = lga( pp*pm/(l*p) );//-K2/(4.*p*p);
 res = ( fabs(p/k0)>1e5 ) ? (K2/SQR(2.*p))*(  //1. // dangerous!
                                                //+ k0/p
                                                //+ (k*k+7.*k0*k0)/(8.*p*p)
                                                //+ k0*(k*k+3.*k0*k0)/(4.*CUBE(p))
1 + (pow(k,8) + 46*pow(k,6)*
       pow(k0,2) + 
      256*pow(k,4)*pow(k0,4) + 
      466*pow(k,2)*pow(k0,6) + 
      511*pow(k0,8))/(1280.*pow(p,8)) + 
   (3*pow(k,6)*k0 + 
      31*pow(k,4)*pow(k0,3) + 
      73*pow(k,2)*pow(k0,5) + 
      85*pow(k0,7))/(192.*pow(p,7)) + 
   (pow(k,6) + 29*pow(k,4)*pow(k0,2) + 
      99*pow(k,2)*pow(k0,4) + 
      127*pow(k0,6))/(256.*pow(p,6)) + 
   (k0*pow(pow(k,2) + 3*pow(k0,2),2))/
    (16.*pow(p,5)) + 
   (pow(k,4) + 16*pow(k,2)*pow(k0,2) + 
      31*pow(k0,4))/(48.*pow(p,4)) + 
   (k0*(pow(k,2) + 3*pow(k0,2)))/
    (4.*pow(p,3)) + 
   (pow(k,2) + 7*pow(k0,2))/
    (8.*pow(p,2)) + k0/p
                                               ) : res;//*/
  return res;
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
      pm = (fabs(pm)<1e-1*fabs(km-k)) ? pd : pm;     // p=km (!)
      if (k0>3.*k) {
        temp +=  lga(pp/pm)*(                         // refl. log
                                      F_123(p,q,r)      // 1A[p,q]
                                    + F_123(q,p,r)      // 1A[q,p]
                                    + F_345(r,p,q)
                                    + F_345(r,q,p)
                                    );
      };
      return .5*temp/r;
  },  kp-p,km  )(y); },  k, km    )(x);              //*/
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
                                 F_123(p,q,r)         // 1B[p,q]
                               + F_123(q,p,r)         // 1B[q,p]
                               + F_345(r,p,q)
                               + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  km-p,min(p,kp-p) )(y); }, .5*km,km   )(x);      //*/
  //\
  return _1B;
  _1D =                                               // corner, +1C
  remap([&](double pq, double pqd) {                  // p-q=[k-min(k,km),k+min(k,km)]
    return remap([&](double r, double rd) {           // r=[|pq-k|,min(k,km)]
      double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
      pm = km-p; qm = km-q;
      pp = kp-p; qp = kp-q;
      double l=k0-p;
      double v=k0-q;
      //\
      Padding for log singularities:
      double h = fabs( min(k,km)-fabs(pq-k) );        // p=kp,q=km (!!)
      //double h = fabs( 2.*min(k,km) );        // p=kp,q=km (!!)
      if (pq>k) { pp = (fabs(pp)<1e-1*h) ? rd/2. : pp; }
      if (pq<k) { qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; }
      double lg1 = lga( pp*qp/(p*q) ),
             lg2 = lga( pp*v/(qm*p) ),
             lg3 = lga( p*q/( l*v ) );

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lg1
                                    + (.5+f(r,s3))*lg2
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 1D[p,q]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lg1
                                    + (.5+f(r,s3))*lg2
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 1C[q,p]
      temp += pow(v/k0,n)*F_14(l,p)*( (.5+f(q,s5))*lg1
                                    + (.5+f(r,s3))*lg2
                                    ); // from (*)
      temp += pow(v/k0,m)*F_25(l,p)*( (.5+f(q,s4))*lg1
                                    + (.5+f(r,s3))*lg2
                                    ); //      (**)
      // Adjust logs:
      lg1 += lg3; lg2 += lg3;

      temp -= pow(v/k0,n)*F_14(l,p)*( (.5+f(v,s2))*lg1
                                    - (.5+f(r,s3))*lg2
      // (***)                      - (.5+f(q,s5))*lga( pm*qm/(p*q) )
                                    ); // 4C[l,v]
      temp -= pow(v/k0,m)*F_25(l,p)*( (.5+f(v,s1))*lg1
                                    - (.5+f(r,s3))*lg2
      // (****)                     - (.5+f(q,s4))*lga( pm*qm/(p*q) )
                                    ); // 4B[v,l]
      temp -= pow(q/k0,n)*F_14(p,l)*( (.5+f(v,s5))*lg1
                                    - (.5+f(r,s3))*lg2
                                    ); //      (***)
      temp -= pow(q/k0,m)*F_25(p,l)*( (.5+f(v,s4))*lg1
                                    - (.5+f(r,s3))*lg2
                                    ); //      (****)
      return .25*temp/r;
  },  fabs(pq-k),min(k,km)    )(y);
  }, k-min(k,km),k+min(k,km)  )(x);                   //*/
  //\
  return _1D;

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
  remap([&](double p, double pd) {                    // p=[kp,+inf)
    pm = km-p; pp = kp-p;
    pp = (fabs(pp)<k0) ? pd : pp;                     // p=kp (!)
    return remap([&](double q, double qd) {           // q=(-inf,km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      //double Lp = lga(pm/pp), Lq = lga(qm/qp);
      //if (+p/k>1e4) { Lp = (k/q)*(1.+.5*k0/p+(k0*k0+k*k/3.)/(4*p*p)); }
      //if (-q/k>1e4) { Lq = (k/q)*(1.+.5*k0/q+(k0*k0+k*k/3.)/(4*q*p)); }
      temp += lga(pm*qm/(pp*qp))*(
                                   F_123(p,q,r)       // 2A[p,q]
                                 + F_123(q,p,r)       // 3A[q,p]
                                 + F_345(r,p,q)
                                 + F_345(r,q,p)
                                 );
      return .5*temp/r;
  },  km-p  )(y); },  kp )(x);                        //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) {                    // p=[kp,+inf)
    pp = kp-p;
    pp = (fabs(pp)<k0) ? pd : pp;                     // p=kp (!)
    return remap([&](double q, double qd) {           // q=[km-p,kp-p]
      double temp = 0., r = k0-p-q;
      qp = kp-q;
      temp += lga(p*q/(pp*qp))*(
                                 F_123(p,q,r)         // 2B[p,q]
                               + F_123(q,p,r)         // 3B[q,p]
                               + F_345(r,p,q)
                               + F_345(r,q,p)
                               );
      return .5*temp/r;
  },  km-p,kp-p  )(y); },  kp )(x);                   //*/
  //\
  return _2B;
  _2C =
  remap([&](double p, double pd) {                    // p=[km,kp]
    pm = km-p; pp = kp-p;
    pm = (fabs(pm)<1e-1*k) ? pd : pm;                 // p=km (!)
    pp = (fabs(pp)<1e-1*k) ? pd : pp;                 // p=kp (!)
    return remap([&](double q, double qd) {           // q=(-inf,km-p]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;

       temp += pow(q/k0,n)*F_14(p,l)*( (-f(fabs(q),s2))*lga( pm*qm/(p*q) )
                                     + (+f(fabs(r),s3))*lga( pm*v/(qp*p) )
                                     - .5*L_div(q,v,qp,qm)
                                     + ( (n==2) ? K2/(8.*q*q) : 0.)
       // (*)                        - (.5+f(v,s5))*lga( pp*qp/(l*v) )
                                     ); // 2C[p,q]
       temp += pow(q/k0,m)*F_25(p,l)*( (-f(fabs(q),s1))*lga( pm*qm/(p*q) )
                                     + (+f(fabs(r),s3))*lga( pm*v/(qp*p) )
                                     - .5*L_div(q,v,qp,qm)
                                     + ( (m==2) ? K2/(8.*q*q) : 0.)
       // (**)                       - (.5+f(v,s4))*lga( pp*qp/(l*v) )
                                     ); // 3C[q,p]
       temp += pow(v/k0,n)*F_14(l,p)*( (-f(fabs(q),s5))*lga( pm*qm/(p*q) )
                                     + (+f(fabs(r),s3))*lga( pm*v/(qp*p) )
                                     - .5*L_div(v,q,qm,qp)
                                     + ( (n==2) ? K2/(8.*v*v) : 0.)
                                     ); // from (*)
       temp += pow(v/k0,m)*F_25(l,p)*( (-f(fabs(q),s4))*lga( pm*qm/(p*q) )
                                     + (+f(fabs(r),s3))*lga( pm*v/(qp*p) )
                                     - .5*L_div(v,q,qm,qp)
                                     + ( (m==2) ? K2/(8.*v*v) : 0.)
                                     ); //      (**)
       return .5*temp/r;
  },  km-p  )(y); },  km,kp )(x);                     //*/
  //\
  return _2C;
  _2D =
  remap([&](double p, double pd) {                    // p=[km,kp]
    return remap([&](double q, double qd) {           // q=[km-p,0]
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      //double q_=q;
      double q_ = (fabs(q)<1e-1*fabs(km-p)) ? qd : q; // q=0 (!)
      double l=k0-p, p2=p*p;
      double v=k0-q, q2=q*q;

      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r,s3))*lga( q_*v/( qp*qm ) )
                                        ); // 2D[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r,s3))*lga( q_*v/( qp*qm ) ) 
                                        ); // 3D[p,q]

      temp += pow(q/k0,m)*F_25(p,l)*( (m==2) ? K2/(8.*q*q) : 0.); // Case: m=2
      temp += pow(v/k0,m)*F_25(l,p)*( (m==2) ? K2/(8.*v*v) : 0.);
      return .5*temp/r;
  },  km-p, 0.  )(y); },  km,kp    )(x);              //*/
  //\
  return _2D;

  //  #] from '2','3'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) {                    // p=[kp,+inf)
    pm = km-p; pp = kp-p;
    pp = (fabs(pp)<k0) ? pd : pp;                     // p=kp (!)
    return remap([&](double q, double qd) {           // q=[p,+inf)
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      //qp = (fabs(qp)<k0) ? qd : qp; // q=kp (!)
      //if (p>q) { pp = (fabs(pp)<1e-1*h) ? pqd/2. : pp; }
      //if (p<q) { qp = (fabs(qp)<1e-1*h) ? pqd/2. : qm; }
      temp += lga(pm/(pp))*(
                                   F_123(p,q,r)       // 4A[p,q]
                                 + F_123(q,p,r)       // 4A[q,p]
                                 + F_345(r,p,q)
                                 + F_345(r,q,p)
                                 );
      return .5*temp/r;
  },  kp  )(y); },  kp    )(x);                        //*/
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
      l = (fabs(l)<1e-1*k0) ? pd : l;
      /*double LL = lga( pp*pm/(l*p) );//-K2/(4.*p*p);
      LL = ( fabs(p/k0)>1e4 ) ? (K2/SQR(2.*p))*(
                                              1.+ k0/p
                                                + (k*k+7.*k0*k0)/(8.*p*p)
                                                + k0*(k*k+3.*k0*k0)/(4.*CUBE(p))
                                               ) : LL;*/

      temp += pow(p/k0,n)*F_14(q,v)*( (+f(fabs(p),s2))*lga( pm*qm/(p*q) )
                                    + (-f(fabs(r),s3))*lga( qm*l/(pp*q) )
                                    + .5*L_div(p,l,pp,pm)
                                    - ( (n==2) ? K2/(8.*p*p) : 0.)
      // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                    ); // 4C[q,p]
      temp += pow(p/k0,m)*F_25(q,v)*( (+f(fabs(p),s1))*lga( pm*qm/(p*q) )
                                    + (-f(fabs(r),s3))*lga( qm*l/(pp*q) )
                                    + .5*L_div(p,l,pp,pm)
                                    //+ (.5)*lga( pm*pp/(p*l) )
      // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                    - ( (m==2) ? K2/(8.*p*p) : 0.)
                                    ); // 4B[p,q]
      temp += pow(l/k0,n)*F_14(v,q)*( (+f(fabs(p),s5))*lga( pm*qm/(p*q) )
                                    + (-f(fabs(r),s3))*lga( qm*l/(pp*q) )
                                    + .5*L_div(l,p,pm,pp)
                                    - ( (n==2) ? K2/(8.*l*l) : 0.)
                                    ); // from (*)
      temp += pow(l/k0,m)*F_25(v,q)*( (+f(fabs(p),s4))*lga( pm*qm/(p*q) )
                                    + (-f(fabs(r),s3))*lga( qm*l/(pp*q) )
                                    + .5*L_div(l,p,pm,pp)
                                    - ( (m==2) ? K2/(8.*l*l) : 0.)
                                    ); //      (**)
      return .5*temp/r;
  },  k0  )(y); },  km,kp    )(x);                    //*/
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
        //qm = (fabs(qm)<1e-1*k) ? qd : qm;               // q=km (!)
        double l=k0-p; l = (fabs(l)<1e-1*fabs(k)) ? pd : l;
        double v=k0-q;
        temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
        // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                      ); // 4C[q,p]
        temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
        // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                      ); // 4B[p,q]
        temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                      ); // from (*)
        temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                      ); //      (**)

        temp += ( pow(l/k0,n)*F_14(v,q)
                + pow(p/k0,n)*F_14(q,v) )*( (.5+f(r,s3))*lga( l*p/( pp*pm ) )
                                          ); // 2D[v,l]
        temp += ( pow(l/k0,m)*F_25(v,q)
                + pow(p/k0,m)*F_25(v,q) )*( (.5+f(r,s3))*lga( l*p/( pp*pm ) )
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
        if (p>2.*kp-km) { qm = (fabs(qm)<1e-1*h) ? qd : qm; }              // q=km (!)
        double l = k0-p, v=k0-q; 
        l  = ( fabs(l)<1e-1*km) ? pd : l;
        pp = (fabs(pp)<1e-1*km) ? pd : pp;
        temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
        // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                      ); // 4C[q,p]
        temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
        // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                      ); // 4B[p,q]
        temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                      ); // from (*)
        temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                      ); //      (**)

        temp -= pow(l/k0,n)*F_14(v,q)*( (.5+f(l,s2))*lga( pm*qm/(v*l) )
                                      + (.5+f(r,s3))*lga( pp*q/(qm*p) )
        // (***)                      - (.5+f(p,s5))*lga( pp*qp/(q*p) )
                                      ); // 1D[v,l]
        temp -= pow(l/k0,m)*F_25(v,q)*( (.5+f(l,s1))*lga( pm*qm/(v*l) )
                                      + (.5+f(r,s3))*lga( pp*q/(qm*p) )
        // (****)                     - (.5+f(p,s4))*lga( pp*qp/(q*p) )
                                      ); // 1C[l,v]
        temp -= pow(p/k0,n)*F_14(q,v)*( (.5+f(l,s5))*lga( pm*qm/(v*l) )
                                      + (.5+f(r,s3))*lga( pp*q/(qm*p) )
                                      ); //      (***)
        temp -= pow(p/k0,m)*F_25(q,v)*( (.5+f(l,s4))*lga( pm*qm/(v*l) )
                                      + (.5+f(r,s3))*lga( pp*q/(qm*p) )
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
        temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
        // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                      ); // 4C[q,p]
        temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
        // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                      ); // 4B[p,q]
        temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                      ); // from (*)
        temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                      + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                      ); //      (**)

        temp += ( pow(l/k0,n)*F_14(v,q)
                + pow(p/k0,n)*F_14(q,v) )*( (.5+f(r,s3))*lga( l*p/( pp*pm ) )
                                          ); // 2D[v,l]
        temp += ( pow(l/k0,m)*F_25(v,q)
                + pow(p/k0,m)*F_25(v,q) )*( (.5+f(r,s3))*lga( l*p/( pp*pm ) )
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

      //double h = fabs( k-fabs(r) );
      //if (r>0) { qm = (fabs(qm)<1e-1*h) ? pqd/2. : qm; }           // q=km (!)
      //if (r<0) { pp = (fabs(pp)<1e-1*h) ? pqd/2. : pp; }           // q=km (!)
      double h = fabs( k-(pq) );
      qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lga( qm*qp/q2 )
                                    - (.5+f(v,s5))*lga( qm*qp/v2 ) ); // 4D[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*pp/p2 )
                                    - (.5+f(l,s4))*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*pp/p2 )
                                    - (.5+f(l,s5))*lga( pm*pp/l2 ) ); // 4D[q,p]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lga( qm*qp/q2 )
                                    - (.5+f(v,s4))*lga( qm*qp/v2 ) );

      temp -= pow(v/k0,n)*F_14(l,p)*( (.5+f(v,s2))*lga( qm*qp/v2 )
                                    - (.5+f(q,s5))*lga( qm*qp/q2 ) ); // 4D[l,v]
      temp -= pow(l/k0,m)*F_25(v,q)*( (.5+f(l,s1))*lga( pm*pp/l2 )
                                    - (.5+f(p,s4))*lga( pm*pp/p2 ) );

      temp -= pow(l/k0,n)*F_14(v,q)*( (.5+f(l,s2))*lga( pm*pp/l2 )
                                    - (.5+f(p,s5))*lga( pm*pp/p2 ) ); // 4D[v,l]
      temp -= pow(v/k0,m)*F_25(l,p)*( (.5+f(v,s1))*lga( qm*qp/v2 )
                                    - (.5+f(q,s4))*lga( qm*qp/q2 ) );
      return .25*temp/r;
  //},  0.,k-fabs(r) )(y); }, -k,k    )(x);                  //*/
  },  0.,k-(pq) )(y); }, 0.,k    )(x);                  //*/
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

  res += _1A + _1B + _1D + _2A + _2B + _2C + _2D + _4A + _4B + _4D;

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
              + pow(l/k0,n)*F_14(v,q) )*( (.5+f(r,s3))*lg1
                                       ); // 2A[p,q]
      temp += ( pow(p/k0,m)*F_25(q,v)
              + pow(l/k0,m)*F_25(v,q) )*( (.5+f(r,s3))*lg1
                                        ); // 3A[q,p]
      return .5*temp/r;
  },  km  )(y); },  km,0.    )(x);
  //
  _2A +=
  remap([&](double p, double pd) {                    // p=[0,kp]
    return remap([&](double q, double qd) {           // q=(-inf,km]
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double lg1 = lga( p*l/(pp*pm) );
      temp += ( pow(p/k0,n)*F_14(q,v)
              + pow(l/k0,n)*F_14(v,q) )*( (.5+f(r,s3))*lg1
                                        ); // 2A[p,q]
      temp += ( pow(p/k0,m)*F_25(q,v)
              + pow(l/k0,m)*F_25(v,q) )*( (.5+f(r,s3))*lg1
                                        ); // 3A[p,q]
      return .5*temp/r;
  },  km-p  )(y); },  0.,kp    )(x);                  //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) {                    // p=[0,k0]
    return remap([&](double q, double qd) {           // q=[km,+inf)
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double h = fabs( p );
      qm = (fabs(qm)<1e-1*h) ? qd : qm;

      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 2B[p,q]
      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 3B[q,p]
      temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                    ); // from (*)
      temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                    ); //      (**)
      return .5*temp/r;
  },  km-p,km  )(y); },  0.,k0    )(x);
  //
  _2B +=
  remap([&](double p, double pd) {                    // p=[0,k0]
    return remap([&](double q, double qd) {           // q=[km,+inf)
      double temp = 0., r = k0-p-q;
      qm = km-q; qp = kp-q;
      pm = km-p; pp = kp-p;
      double l=k0-p, v=k0-q;
      double h = fabs( kp-k0 );
      pp = (fabs(pp)<1e-1*h) ? pd : pp;

      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 2B[p,q]
      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 3B[q,p]
      temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                    ); // from (*)
      temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
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
      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
      // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                    ); // 2B[p,q]
      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
      // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                    ); // 3B[q,p]
      temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                    ); // from (*)
      temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                    + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                    ); //      (**)

      temp -= pow(l/k0,n)*F_14(v,q)*( (.5+f(l,s2))*lga( pm*qm/(l*v) )
                                    - (.5+f(r,s3))*lga( qm*p/(pp*v) )
                                    ); // 4B'[v,l]
      temp -= pow(l/k0,m)*F_25(v,q)*( (.5+f(l,s1))*lga( pm*qm/(l*v) )
                                    - (.5+f(r,s3))*lga( qm*p/(pp*v) )
                                    ); // 4C'[l,v]
      temp -= pow(p/k0,n)*F_14(q,v)*( (.5+f(l,s5))*lga( pm*qm/(l*v) )
                                    - (.5+f(r,s3))*lga( qm*p/(pp*v) )
                                    );
      temp -= pow(p/k0,m)*F_25(q,v)*( (.5+f(l,s4))*lga( pm*qm/(l*v) )
                                    - (.5+f(r,s3))*lga( qm*p/(pp*v) )
                                    );
      return .25*temp/r;
  },  fabs(pq-k),-km  )(y); }, k+km,k-km  )(x);       //*/
  //\
  return _2B;

  //  #] from '2'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) {                    // p=[kp,+inf)
    return remap([&](double q, double qd) {           // q=[kp,+inf)
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

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lga( qm*qp/q2 )
                                    - (.5+f(v,s5))*lga( qm*qp/v2 ) ); // 6B[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*pp/p2 )
                                    - (.5+f(l,s4))*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*pp/p2 )
                                    - (.5+f(l,s5))*lga( pm*pp/l2 ) ); // 7B[q,p]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lga( qm*qp/q2 )
                                    - (.5+f(v,s4))*lga( qm*qp/v2 ) );
      temp += lga(qm*pm/( p*q ))*(
                                  F_123(p,q,r)
                                + F_123(q,p,r)
                                + F_345(r,p,q)
                                + F_345(r,q,p)
                               ); 
      return .5*temp/r;
  },  kp,p  )(y); },  kp    )(x);                     //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) {                    // q=[kp,k0]
    return remap([&](double p, double pd) {           // p=[kp,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r,s3))*lga( q*v/( qp*qm ) )
                                            ); // 4B[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r,s3))*lga( q*v/( qp*qm ) ) 
                                            ); // 4C[q,p]
      return .5*temp/r;
  },  kp  )(y); },  k0+km,kp    )(x);
  //
  _4B +=
  remap([&](double q, double qd) {                    // q=[0,k0]
    return remap([&](double p, double pd) {           // p=[kp,+inf)
      double temp = 0., r = k0-p-q;
      pm = km-p; pp = kp-p;
      qm = km-q; qp = kp-q;
      double l=k0-p, v=k0-q;
      temp += ( pow(q/k0,n)*F_14(p,l)
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r,s3))*lga( q*v/( qp*qm ) )
                                        ); // 4B[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r,s3))*lga( q*v/( qp*qm ) ) 
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
              + pow(v/k0,n)*F_14(l,p) )*( (.5+f(r,s3))*lga( q*v/( qp*qm ) )
                                        ); // 4B[p,q]
      temp += ( pow(q/k0,m)*F_25(p,l)
              + pow(v/k0,m)*F_25(l,p) )*( (.5+f(r,s3))*lga( q*v/( qp*qm ) ) 
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
    return remap([&](double p, double pd) {           // p=(-inf,km]
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

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lga( qm*qp/q2 )
                                    - (.5+f(v,s5))*lga( qm*qp/v2 ) ); // 6B[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*pp/p2 )
                                    - (.5+f(l,s4))*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*pp/p2 )
                                    - (.5+f(l,s5))*lga( pm*pp/l2 ) ); // 7B[q,p]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lga( qm*qp/q2 )
                                    - (.5+f(v,s4))*lga( qm*qp/v2 ) );

      temp += lga(qp*pp/(p*q))*(
                                  F_123(p,q,r)
                                + F_123(q,p,r)
                                + F_345(r,p,q)
                                + F_345(r,q,p)
                               );

      return .5*temp/r;
  },  q,km  )(y); },  km    )(x);                     //*/
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

      temp += lga( pm*qm/(p*q) )*( F_123(p,q,r) + F_345(r,p,q) ); // 5B[p,q]

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

      temp += lga(qp*pp/( p*q ))*(
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
  remap([&](double pq, double pqd) {                  // p-q=[k,k-km]
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

      temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lga( qm*qp/q2 )
                                    - (.5+f(v,s5))*lga( qm*qp/v2 ) ); // 6B[p,q]
      temp += pow(p/k0,m)*F_25(q,v)*( (.5+f(p,s1))*lga( pm*pp/p2 )
                                    - (.5+f(l,s4))*lga( pm*pp/l2 ) );

      temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*pp/p2 )
                                    - (.5+f(l,s5))*lga( pm*pp/l2 ) ); // 7B[q,p]
      temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lga( qm*qp/q2 )
                                    - (.5+f(v,s4))*lga( qm*qp/v2 ) );
      return .25*temp/r;
  },  k-pq,pq-k  )(y); }, k    )(x);                  //*/
  //
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
  _2A + _2B + _4A + _4B + _4D + _5A + _5B + _6A + _6B + _6C;
  }

  // still don't have a good way to cater for NaNs
  //if ( isinf(res)||isnan(res) ) { return -1e6;}//for inspection
  //if ( isinf(res)||isnan(res) ) { return 0.;}
  //else 
  { return .5*(K2)*res/k; }
  //else { return res/k; }
}
