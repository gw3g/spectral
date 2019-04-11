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

  double F_123(double,double,double); // Cut: f1f2f3/f0   \_ real
  double F_345(double,double,double); //      f3f4f5/f0   /
  double F_14(double,double);         //        f1f4/f0   \_ virtual
  double F_25(double,double);         //        f2f5/f0   /

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
    double epsabs = 1e-2, epsrel = 1e-2;
    size_t limit = 1e5;

    quad wsp1(limit);
    quad wsp2(limit);

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return (this->integrand)(x,y);
          } );
      gsl_integration_qag( inner, .0+1e-10,1., epsabs, epsrel, 
                           limit, 6, wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, .0+1e-10,1., epsabs, epsrel, 
                         limit, 6, wsp2, &res, &err  );//*/

    //return ( res -k0*k0/8. )*CUBE(OOFP);
    return ( res )*CUBE(OOFP);
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
  remap([&](double p, double pd) { // [.5kp,km]
       pm = km-p; pp = kp-p;
       pm = (fabs(pm)<1e-1*fabs(km-.5*kp)) ? pd : pm; // p=km (!)
       return remap([&](double q, double qd) { // [kp-p,p]
          double temp = 0., r=k0-p-q;
          qm = km-q; qp = kp-q;
          if (k0>3.*k) {
            temp +=  lga(pp*qp/(pm*qm))*(
                                          F_123(p,q,r) // 1A[p,q]
                                        + F_123(q,p,r) // 1A[q,p]
                                        + F_345(r,p,q)
                                        + F_345(r,q,p)
                                        );
          };
          return .5*temp/r;
  },  kp-p,p  )(y); },  .5*kp, km    )(x);//*/
  //
  _1B =
  remap([&](double p, double pd) { // [.5km,km]
       pm = km-p; 
       pm = (fabs(pm)<1e-1*( .5*km )) ? pd : pm; // p=km (!)
       return remap([&](double q, double qd) { // [km-p,min(p,kp-p)]
          double temp = 0., r = k0-p-q;
          qm = km-q;
          temp += lga(p*q/(pm*qm))*(
                                       F_123(p,q,r) // 1B[p,q]
                                     + F_123(q,p,r) // 1B[q,p]
                                     + F_345(r,p,q)
                                     + F_345(r,q,p)
                                     );
          return .5*temp/r;
  },  km-p,min(p,kp-p) )(y); }, .5*km,km   )(x);//*/
  //
  _1D = // corner, 1D+C
  remap([&](double pq, double pqd) { // [k-min(k,km),k+min(k,km)]
      return remap([&](double r, double rd) { // [|pq-k|,min(k,km)]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          pm = km-p; qm = km-q;
          pp = kp-p; qp = kp-q;
          double l=k0-p;
          double v=k0-q;
          //\
          Padding for log singularities:
          double h = fabs( min(k,km)-fabs(pq-k) );
          if (pq>k) { pp = (fabs(pp)<1e-1*h) ? rd/2. : pp; }
          if (pq<k) { qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; }

          temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lga( pp*qp/(p*q) )
                                        + (.5+f(r,s3))*lga( pp*v/(qm*p) )
          // (*)                        - (.5+f(v,s5))*lga( pm*qm/(l*v) )
                                        ); // 1D[p,q]
          temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lga( pp*qp/(p*q) )
                                        + (.5+f(r,s3))*lga( pp*v/(qm*p) )
          // (**)                       - (.5+f(v,s4))*lga( pm*qm/(l*v) )
                                        ); // 1C[q,p]
          temp += pow(v/k0,n)*F_14(l,p)*( (.5+f(q,s5))*lga( pp*qp/(p*q) )
                                        + (.5+f(r,s3))*lga( pp*v/(qm*p) )
                                        ); // from (*)
          temp += pow(v/k0,m)*F_25(l,p)*( (.5+f(q,s4))*lga( pp*qp/(p*q) )
                                        + (.5+f(r,s3))*lga( pp*v/(qm*p) )
                                        ); //      (**)

          temp -= pow(v/k0,n)*F_14(l,p)*( (.5+f(v,s2))*lga( pp*qp/(l*v) )
                                        - (.5+f(r,s3))*lga( pp*q/(qm*l) )
          // (***)                      - (.5+f(q,s5))*lga( pm*qm/(p*q) )
                                        ); // 4C[l,v]
          temp -= pow(v/k0,m)*F_25(l,p)*( (.5+f(v,s1))*lga( pp*qp/(l*v) )
                                        - (.5+f(r,s3))*lga( pp*q/(qm*l) )
          // (****)                     - (.5+f(q,s4))*lga( pm*qm/(p*q) )
                                        ); // 4B[v,l]
          temp -= pow(q/k0,n)*F_14(p,l)*( (.5+f(v,s5))*lga( pp*qp/(l*v) )
                                        - (.5+f(r,s3))*lga( pp*q/(qm*l) )
                                        ); //      (***)
          temp -= pow(q/k0,m)*F_25(p,l)*( (.5+f(v,s4))*lga( pp*qp/(l*v) ) 
                                        - (.5+f(r,s3))*lga( pp*q/(qm*l) ) 
                                        ); //      (****)
          return .25*temp/r;
  },  fabs(pq-k),min(k,km)  )(y); }, k-min(k,km),k+min(k,km)  )(x); //*/\
  return _1D;

  //  #] from '1'\
      #[ from '2' (& by sym, '3')

  _2A =
  remap([&](double p, double pd) { // [kp,+inf)
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // (-inf,km-p]
        double temp = 0., r = k0-p-q;
        qm = km-q; qp = kp-q;
        temp += lga(pm*qm/(pp*qp))*(
                                     F_123(p,q,r) // 2A[p,q]
                                   + F_123(q,p,r) // 2A[q,p]
                                   + F_345(r,p,q)
                                   + F_345(r,q,p)
                                   );
        return .5*temp/r;
  },  km-p  )(y); },  kp )(x); //*/
  //\
  return _2A;
  _2B =
  remap([&](double p, double pd) { // [kp,+inf)
      pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // [km-p,kp-p]
        double temp = 0., r = k0-p-q;
        qp = kp-q;
        temp += lga(p*q/(pp*qp))*(
                                          F_123(p,q,r) // 2B[p,q]
                                        + F_123(q,p,r) // 2B[q,p]
                                        + F_345(r,p,q)
                                        + F_345(r,q,p)
                                 );
        return .5*temp/r;
  },  km-p,kp-p  )(y); },  kp )(x); //*/
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

          temp += pow(q/k0,n)*F_14(p,l)*( (.5+f(q,s2))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( pm*v/(qp*p) )
          // (*)                        - (.5+f(v,s5))*lga( pp*qp/(l*v) )
                                        ); // 2C[p,q]
          temp += pow(q/k0,m)*F_25(p,l)*( (.5+f(q,s1))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( pm*v/(qp*p) )
          // (**)                       - (.5+f(v,s4))*lga( pp*qp/(l*v) )
                                        //+ K2/(8.*q*q) 
                                        ); // 3C[q,p]

          temp += pow(v/k0,n)*F_14(l,p)*( (.5+f(q,s5))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( pm*v/(qp*p) )
                                        ); // from (*)
          temp += pow(v/k0,m)*F_25(l,p)*( (.5+f(q,s4))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( pm*v/(qp*p) )
                                        //+ K2/(8.*v*v) 
                                        ); //      (**)

          return .5*temp/r;
  },  km-p  )(y); },  km,kp )(x); //*/
  //\
  return _2C;
  _2D =
  remap([&](double p, double pd) { // [km,kp]
      return remap([&](double q, double qd) { // [km-p,0]
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
          return .5*temp/r;
  },  km-p, 0.  )(y); },  km,kp    )(x); //*/\
  return _2D;

  //  #] from '2','3'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) { // [kp,+inf)
      pm = km-p; pp = kp-p;
      pp = (fabs(pp)<k0) ? pd : pp; // p=kp (!)
      return remap([&](double q, double qd) { // [p,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          //qp = (fabs(qp)<k0) ? qd : qp; // q=kp (!)
          temp += lga(pm*qm/(pp*qp))*(
                                       F_123(p,q,r) // 4A[p,q]
                                     + F_123(q,p,r) // 4A[q,p]
                                     + F_345(r,p,q)
                                     + F_345(r,q,p)
                                     );
          return .5*temp/r;
  },  p  )(y); },  kp    )(x); //*/
  //\
  return _4A;
  _4B =
  remap([&](double q, double qd) { // [km,kp]
      qm = km-q; qp = kp-q;
      qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=km (!)
      qp = (fabs(qp)<1e-1*k) ? qd : qp; //  =kp (!)
      return remap([&](double p, double pd) { // [k0,+inf)
          double temp = 0., r = k0-p-q;
          pm = km-p; pp = kp-p;
          double l = k0-p, v=k0-q;
          double LL = lga( pp*pm/(l*p) );//-K2/(4.*p*p);
          LL = ( p>1e7 ) ? (K2/SQR(2.*p))*( 1.+ k0/p + (k*k+7.*k0*k0)/(8.*p*p) 
              + k0*(k*k+3.*k0*k0)/(4.*CUBE(p)) )
          : LL;
          temp += pow(p/k0,n)*F_14(q,v)*( (.5+f(p,s2))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( qm*l/(pp*q) )
          // (*)                        - (.5+f(l,s5))*lga( pp*qp/(l*v) )
                                        ); // 4C[q,p]
          temp += pow(p/k0,m)*F_25(q,v)*( 
                                          f(fabs(p),s1)*lga( pm*qm/(p*q) )
                                        - f(fabs(r),s3)*lga( qm*l/(pp*q) )
                                        + .5*LL
                                        //+ (.5)*lga( pm*pp/(p*l) )
          // (**)                       - (.5+f(l,s4))*lga( pp*qp/(l*v) )
                                        //- K2/(8.*p*p) 
                                        ); // 4B[p,q]
          temp += pow(l/k0,n)*F_14(v,q)*( (.5+f(p,s5))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( qm*l/(pp*q) )
                                        ); // from (*)
          temp += pow(l/k0,m)*F_25(v,q)*( (.5+f(p,s4))*lga( pm*qm/(p*q) )
                                        + (.5+f(r,s3))*lga( qm*l/(pp*q) )
          //                              - K2/(8.*l*l) 
                                        ); //      (**)
          return .5*temp/r;
  },  k0  )(y); },  km,kp    )(x);
  //\
  return _4B;
  if (k0>3.*k) {
  //
    _4B+=
    remap([&](double q, double qd) { // [km,kp]
        qm = km-q;
        qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=km (!)
        return remap([&](double p, double pd) { // [k0+km-q,k0]
            double temp = 0., r = k0-p-q;
            pm = km-p; pp = kp-p;
            double l=k0-p; l = (fabs(l)<1e-1*fabs(q-km)) ? pd : l;
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
    },  k0+km-q,k0  )(y); },  km,kp    )(x);
  //
    _4B+=
    remap([&](double q, double qd) { // [km,kp]
        qm = km-q;
        qm = (fabs(qm)<1e-1*k) ? qd : qm; // q=km (!)
        return remap([&](double p, double pd) { // [2kp-q,k0+km-q]
            double temp = 0., r = k0-p-q;
            pm = km-p; pp = kp-p;
            double l = k0-p, v=k0-q; //l = (fabs(l)<1e-1*(km-k)) ? pd : l;
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
    },  2.*kp-q,k0+km-q  )(y); },  km,kp    )(x);
  //
  } else if (k0<3.*k) {
  //
    _4B+=
    remap([&](double p, double pd) { // [kp,k0]
        pm = km-p; pp = kp-p;
        pp = (fabs(pp)<1e-1*fabs(k0-kp)) ? pd : pp; // p=kp (!)
        double l = k0-p; l = (fabs(l)<1e-1*fabs(k0-kp)) ? pd : l; // p=k0 (!)
        return remap([&](double q, double qd) { // [k0+km-p,kp]
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
    },  k0+km-p, kp  )(y); },  kp,k0    )(x); //*/
  }
  //\
  return _4B;
  _4D =
  remap([&](double pq, double pqd) { // [max(0,k-km),k]
      return remap([&](double r, double rd) { // [p-q-k,k-p+q]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          pm = km-p; qm = km-q;
          pp = kp-p; qp = kp-q;
          double l=k0-p, p2=p*p, l2=l*l;
          double v=k0-q, q2=q*q, v2=v*v;

          double h = fabs( k-pq );
          qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; // q=km (!)

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
  },  0.,k-pq )(y); }, 0.,k    )(x); //*/
  //
  if (k0<3.*k) { 
  //
    _4D +=
    remap([&](double pq, double pqd) { // [0,k-km]
        return remap([&](double r, double rd) { // [km,k-pq]
            double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
            pm = km-p; qm = km-q;
            pp = kp-p; qp = kp-q;
            double h = fabs(k-pq-km);
            double l=k0-p, p2=p*p, l2=l*l;
            double v=k0-q, q2=q*q, v2=v*v;
            qm = (fabs(qm)<1e-1*h) ? rd/2. : qm; // q=km (!)

            temp -= lga(qp*pp/(p*q))*(
                                       F_123(p,q,r) 
                                     + F_123(q,p,r) 
                                     + F_345(r,p,q)
                                     + F_345(r,q,p)
                                     );
            return .25*temp/r;
    },  km,k-pq  )(y); }, 0.,k-km    )(x);
  //
  }
  //\
  return _4D;

  res += _1A + _1B + _1D + _2A + _2B + _2C + _2D + _4A + _4B + _4D;

  } else {

  double _2A,
         _4A, _4B,      _4D,
         _5A, _5B,
         _6A, _6B, _6C; // 7 by reflection

  //  #[ from '2'

  _2A =
  remap([&](double p, double pd) { // [km,0]
      pm = km-p; pp = kp-p;
      return remap([&](double q, double qd) { // [km,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  km  )(y); },  km,0.    )(x);

  _2A +=
  remap([&](double p, double pd) { // [0,k0]
      pm = km-p, pp = kp-p;
      return remap([&](double q, double qd) { // [km,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  km  )(y); },  0.,k0    )(x);

  _2A +=
  remap([&](double p, double pd) { // [k0,kp]
      pm = km-p, pp = kp-p;
      return remap([&](double q, double qd) { // [k0+km-p,+inf)
      double temp = 0., r = k0-p-q;
      temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
      temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
      return temp/r;
  },  k0+km-p  )(y); },  k0,kp    )(x); //*/

  //  #] from '2'\
      #[ from '4'

  _4A =
  remap([&](double p, double pd) { // [kp,+inf)
      return remap([&](double q, double qd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          temp += W_vi(p,q)*( F_123(p,q,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp  )(y); },  kp    )(x); //*/

  _4B =
  remap([&](double q, double qd) { // [kp,k0]
      qm = km-q; qp = kp-q;
      return remap([&](double p, double pd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp  )(y); },  k0,kp    )(x); //*/

  _4B +=
  remap([&](double q, double qd) { // [0,k0]
      qm = km-q; qp = kp-q;
      return remap([&](double p, double pd) { // [kp,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp  )(y); },  0.,k0    )(x);

  _4B +=
  remap([&](double q, double qd) { // [km,0]
      qm = km-q; qp = kp-q;
      return remap([&](double p, double pd) { // [kp-q,+inf)
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  kp-q  )(y); }, km, 0.  )(x); //*/

  _4D =
  remap([&](double q, double qd) { // [0,kp]
      return remap([&](double p, double pd) { // [kp-q,kp]
          double temp = 0., r = k0-p-q;
          qm = km-q; qp = kp-q;
          temp += W_vi(p,q)*( F_123(p,q,r) );
          return temp/r;
  },  kp-q,kp  )(y); },  0.,kp    )(x); //*/

  //  #] from '4'\
      #[ from '5'

  _5A =
  remap([&](double q, double qd) { // (-inf,km]
      qm = km-q; qp = kp-q;
      return remap([&](double p, double pd) { // (-inf,km]
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          return temp/r;
  },  km  )(y); },  km    )(x); //*/

  _5B =
  remap([&](double q, double qd) { // [km,0]
      qm = km-q, qp = kp-q;
      return remap([&](double p, double pd) { // [km,km-q]
          double temp = 0., r = k0-p-q;
          temp += W_vi(p,q)*( F_123(p,q,r) );
          return temp/r;
  },  km,km-q  )(y); },  km,0    )(x); //*/

  //  #] from '5'\
      #[ from '6'

  _6A =
  remap([&](double p, double pd) { // [k,+inf)
      pm = km-p, pp = kp-p;
      return remap([&](double q, double qd) { // [kp-p,km]
          double temp = 0., r = k0-p-q;
          qm = km-q, qp = kp-q;
          temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );
          return temp/r;
  },  kp-p,km  )(y); },  k    )(x); //*/

  _6B =
  remap([&](double pq, double pqd) { // [k-km,+inf)
      return remap([&](double r, double rd) { // [km,-km]
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

          temp -= F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qp/(q*q) );
          temp -= F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pp/(p*p) );
          temp -= F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qp/(q*q) );
          temp -= F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pp/(p*p) );

          return .5*temp/r;
  },  km,-km  )(y); }, k-km    )(x); //*/

  _6B +=
  remap([&](double p, double pd) { // [kp,+inf)
      pm = km-p, pp = kp-p;
      return remap([&](double q, double qd) { // [km-p,km-p+k0]
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
  remap([&](double pq, double pqd) { // [k,k-km]
      return remap([&](double r, double rd) { // [k-pq,pq-k]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          qm = km-q, qp = kp-q;
          pm = km-p, pp = kp-p;
          double h = 2.*fabs( pq-k );

          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
          //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
          //temp += W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          //temp += lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          //temp += lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );

          temp -= F_14(p,k0-p)*(.5+f(q,s2))*lga( qm*qp/(q*q) );
          temp -= F_14(q,k0-q)*(.5+f(p,s2))*lga( pm*pp/(p*p) );
          temp -= F_25(p,k0-p)*(.5+f(q,s1))*lga( qm*qp/(q*q) );
          temp -= F_25(q,k0-q)*(.5+f(p,s1))*lga( pm*pp/(p*p) );
          return .5*temp/r;
  },  k-pq,pq-k  )(y); }, k,k-km    )(x); //*/


  _6B +=
  remap([&](double pq, double pqd) { // [k+km,k-km]
      return remap([&](double r, double rd) { // [|p-q-k|,-km]
          double temp = 0., p=(k0+pq-r)*.5, q=(k0-pq-r)*.5;
          qm = km-q, qp = kp-q;
          pm = km-p, pp = kp-p;
          double h = fabs( km+fabs(pq-k) );

          //pp = (fabs(pp)<1e-1*h) ? rd/2. : pp;
          //qm = (fabs(qm)<1e-1*h) ? rd/2. : qm;
          //pm = (fabs(pm)<1e-1*h) ? rd/2. : pm;
          //qp = (fabs(qp)<1e-1*h) ? rd/2. : qp;
          temp -= W_vi(p,q)*( F_123(p,q,r) + F_123(q,p,r) );
          temp += W_vi(k0-q,k0-p)*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );

          //temp += lga( pp*qp/(q*p) )*( F_123(p,q,r) + F_123(q,p,r) );
          //temp -= lga( pp*qp/((k0-p)*(k0-q)))*( F_123(k0-q,k0-p,-r) + F_123(k0-p,k0-q,-r) );
          temp -= lga( q*(k0-q)/(qp*qm) )*(.5+f(r,s3))*( F_14(p,k0-p) + F_25(p,k0-p) );
          temp -= lga( p*(k0-p)/(pp*pm) )*(.5+f(r,s3))*( F_14(q,k0-q) + F_25(q,k0-q) );

          return .5*temp/r;
  },  fabs(pq-k),-km  )(y); }, k+km,k-km  )(x); //*/

  _6C =
  remap([&](double q, double qd) { // (-inf,-k]
      qm = km-q, qp = kp-q;
      return remap([&](double p, double pd) { // [kp,km-q]
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
