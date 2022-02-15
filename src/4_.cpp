#include "core.hh"
#include "quad.hh"
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

  double F_123(double,double,double); // Cut: f1f2f3/f0
  double eval();

  rho11100(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=4; }
};
// function for MAIN
Master* _11100(int m, int n, int s[3]) {
  Master *R =  new rho11100(m,n,s); return R;
}

/*--------------------------------------------------------------------*/
// thermal weight for (single) cut from this topology

double rho11100::F_123(double p, double q, double r) { 
  double res;

  int s1 = (this->s)[1],
      s2 = (this->s)[2],
      s3 = (this->s)[3];

  int _m = this->m;
  int _n = this->n;

  double fp=f(p-MOT1,s1), fq=f(q-MOT2,s2), fr=f(r-MOT3,s3);

  res = ( ((double) s1*s2*s3)*exp(k0)-1. )*fp*fq*fr;
  res*= pow(k0,-_m-_n)*pow(p,_m)*pow(q,_n) ;

  return res;
}

/*--------------------------------------------------------------------*/
// evaluation step, including the OPE

double rho11100::eval() {
  double res, err;

  double a1=I(0,(this->s)[1],MOT1), b1=I(2,(this->s)[1],MOT1), // tadpole ints
         a2=I(0,(this->s)[2],MOT2), b2=I(2,(this->s)[2],MOT2),
         a3=I(0,(this->s)[3],MOT3), b3=I(2,(this->s)[3],MOT3);

  if ( m==0 && n==0 ) { // (0)
    (this->OPE).T0 = +K2/8.;
    (this->OPE).T2 = +( a1+a2+a3 )*.25*OOFP;
    (this->OPE).T4 = 0.;
  } else {
    cerr << "Case: (m,n)=("<< m << ","<<n<<") out of bounds!\n";
    return 0.;
  }
  //return (this->OPE)();

  double epsabs = 1e-3, epsrel = 0;
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

  return (( res + ( (k0>k) ? -K2/8. : 0. ) ))*CUBE(OOFP);//*/
  // note: subtraction does not introduce any discontinuity
}

/*--------------------------------------------------------------------*/
// integrand manipulations

double rho11100::integrand(double x, double y) 
{
  double res=0.; 

  if (k0>k) {

  res+=
   remap([&](double p, double pd) { return remap([&](double q, double qd) {
      return W_iv(p,q)*F_123(p,q,k0-p-q);
  },  0., k0-p  )(y); },  km,kp    )(x);

  res+=
   remap([&](double p, double pd) { return remap([&](double q, double qd) {
      return W_iv(p,q)*F_123(p,q,k0-p-q);
  },  km-p,   kp     )(y); }, 0.,  km    )(x);

  res+=
   remap([&](double q, double qd) { return remap([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*F_123(p,q,k0-p-q);
      temp += W_iv(q,p)*F_123(q,p,k0-p-q);
      return temp; },  km, kp-q     )(y); }, -k0, .0    )(x); 
  res+=
   remap([&](double q, double qd) { return remap([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*F_123(p,q,k0-p-q);
      temp += W_iv(q,p)*F_123(q,p,k0-p-q);
      return temp;
  },  km, kp-q     )(y); }, -k0    )(x);

  res+= // pq := p-q,
   remap([&](double r, double rd) { return remap([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*F_123(p,q,r);
      temp += W_iv(q,p)*F_123(q,p,r);
      return .5*temp;
  },  0.,k0-2.*km-r )(y); }, -k0,0.    )(x); 
  res+= // pq := p-q,
   remap([&](double r, double rd) { return remap([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*F_123(p,q,r);
      temp += W_iv(q,p)*F_123(q,p,r);
      return .5*temp;
  },  0.,k0-2.*km-r )(y); }, -k0    )(x); 

  } else { // Below the LC:

  res+=
   remap([&](double q, double qd) { return remap([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*F_123(p,q,k0-p-q);
      temp += W_iv(q,p)*F_123(q,p,k0-p-q);
      return temp;
  },  km, kp-q     )(y); }, km    )(x); // 2A,B and 6B,C

  res+=
   remap([&](double p, double pd) { return remap([&](double q, double qd) {
      double temp = 0.;
      temp += W_iv(p,q)*F_123(p,q,k0-p-q);
      temp += W_iv(q,p)*F_123(q,p,k0-p-q);
      return temp;
  },  kp-p, 0.     )(y); }, kp    )(x); // 6,7A

  res+= // pq := p-q,
   remap([&](double r, double rd) { return remap([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*F_123(p,q,r);
      temp += W_iv(q,p)*F_123(q,p,r);
      return .5*temp;
  },  0.,kp )(y); }, km    )(x);  // 4A,D part of B,C??

  res+= // pq := p-q,
   remap([&](double r, double rd) { return remap([&](double pq, double pqd) {
      double temp = 0.;
      double p = .5*(pq-r+k0), q = .5*(-pq-r+k0);
      temp += W_iv(p,q)*F_123(p,q,r);
      temp += W_iv(q,p)*F_123(q,p,r);
      return .5*temp;
  },  km, 0. )(y); }, kp    )(x);  // 5B

  res+=
   remap([&](double p, double pd) { return remap([&](double q, double qd) {
      double temp = 0.;
      temp += W_iv(p,q)*F_123(p,q,k0-p-q);
      temp += W_iv(q,p)*F_123(q,p,k0-p-q);
      return temp;
  },  0.,p-kp     )(y); }, kp    )(x); // part of 4B & c

  res+=
   remap([&](double q, double qd) { return remap([&](double p, double pd) {
      double temp = 0.;
      temp += W_iv(p,q)*F_123(p,q,k0-p-q);
      temp += W_iv(q,p)*F_123(q,p,k0-p-q);
      return temp;
  },  q-km, 0.     )(y); }, km    )(x); // 5A??

  }
  if ( isinf(res)||isnan(res) ) { return 0.;}
  else { return res; }
}
