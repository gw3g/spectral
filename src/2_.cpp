#include "core.hh"
#include "quad.hh"
#include "map.hh"

double F(int nu, int sA, int sB) {
  // here F_m for m={0,1,2}
  double res=0.;
  res += f(kp,sA)*f(km,sB)/km;
  res -= f(km,sA)*f(kp,sB)/kp;
  res *= ((double)sA*sB)*exp(k0)-1.; // = F_0
  if (nu>0) {
    res *= k0;
    res += ( 2.-exp(kp)*((double)(sA+sB)) )*f(kp,sA)*f(kp,sB);
    res += sgn(km)*( f(fabs(km),sA)+f(fabs(km),sB) ); // = F_1
  }
  if (nu>1) {
    res *= k0;
    res -= (km>0)?k:0.;
    res += 2.*(( kp-fabs(km)-lga( (exp(kp)-sB)/(exp(fabs(km))-sB) ) ));
    res += f(fabs(km),sB)*( k0 + f(fabs(km),sA)*exp(fabs(km))*km*((double)(sB-sA)) )*sgn(km);
    res -= f(fabs(kp),sB)*( k0 + f(fabs(kp),sA)*exp(fabs(kp))*kp*((double)(sB-sA)) ); // = F_2
  }

  res *= .0625*OOFP/k; // coeff = 1/16
  return res;
}

/*--------------------------------------------------------------------*/

struct rho11020 : Master {
  double eval();
  rho11020(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=2; }
};
// function(s) for MAIN
Master* _11020(int m, int n, int s[3]) {
  Master *R =  new rho11020(m,n,s); return R;
}
Master* _10120(int m, int n, int s[3]) {
  int s_new[3]; // swap: s_2 <-> s_3
  s_new[0] = s[0];
  s_new[1] = s[1];
  s_new[2] = s[0]*s[1]*s[2];
  Master *R =  new rho11020(m,n,s_new); return R;
}

double rho11020::eval() 
{
  double a1=I(0,(this->s)[1]), b1=I(2,(this->s)[1]), // tadpole ints
         a2=I(0,(this->s)[2]), b2=I(2,(this->s)[2]),
         a4=I(0,(this->s)[4]), b4=I(2,(this->s)[4]);

  if ( m==0 && n==0 ) { // (0)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 =  -a2*.25*OOFP/K2;
    (this->OPE).T4 = 0.;
  } else
  if ( m==1 && n==0 ) { // (1)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 =  -a2*.25*OOFP*k0/K2;
    (this->OPE).T4 = 0.;
  } else
  if ( m==2 && n==0 ) { // (2)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 =  -.75*a2*.25*OOFP*(k0*k0+k*k/3.)/K2;
    (this->OPE).T4 = 0.;
  } else {
    cerr << "Case: (m,n)=("<< m << ","<<n<<") out of bounds!\n";
    return 0.;
  }

  return -F(m,s[1],s[4])*I(n,s[2]); // the minus is from I(..)
}


/*--------------------------------------------------------------------*/

struct rho11010 : Master {
  double integrand(double);           // supported on [0,1]

  double F_14(double,double);         // = f1f4/f0

  double eval();

  rho11010(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=1; }
};
// function for MAIN
Master* _11010(int m, int n, int s[3]) {
  Master *R =  new rho11010(m,n,s); return R;
}
Master* _10110(int m, int n, int s[3]) {
  int s_new[3]; // swap: s_2 <-> s_3
  s_new[0] = s[0];
  s_new[1] = s[1];
  s_new[2] = s[0]*s[1]*s[2];
  Master *R =  new rho11010(m,n,s_new); return R;
}


/*--------------------------------------------------------------------*/
// all the thermal weights for cuts from this topology

double rho11010::F_14(double p, double l) {
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
// evaluation step, including the OPE

double rho11010::eval()
{
  double res, err;

  double a1=I(n,(this->s)[1]), // tadpole ints
         a2=I(n,(this->s)[2]),
         a4=I(n,(this->s)[4]);

  if ( m==0 ) { // (0)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 =  a2*.25*OOFP;
    (this->OPE).T4 = 0.;
  } else
  if ( m==1 ) { // (1)
    (this->OPE).T0 = 0.;
    (this->OPE).T2 =  a2*.125*OOFP*k0;
    (this->OPE).T4 = 0.;
  } else {
    cerr << "Case:  m = "<< m << "  out of bounds!\n";
    return 0.;
  }

  double epsabs = 1e-5, epsrel = 1e-5;
  size_t limit = 1e5;

  quad wsp(limit);

  auto outer = make_gsl_function( [&](double x) {
      return (this->integrand)(x);
  } );
  gsl_integration_qag(  outer, .0+1e-13,1., epsabs, epsrel,
                        limit, 6, wsp, &res, &err  );

  return res*a2*.25*OOFP*pow(k0,m);
}

/*--------------------------------------------------------------------*/
// integrand manipulations

double rho11010::integrand(double x) 
{
  // statistics:
  int s0 = (this->s)[0],
      s1 = (this->s)[1],
      s2 = (this->s)[2],
      s4 = (this->s)[4];

  double res=0.;

  if (k0>k) {
    res +=
    remap([&](double p, double pd) {                    // p=[km,kp]
        return F_14(p,k0-p);
    },  km, kp    )(x);
  } else {
    res +=
    remap([&](double p, double pd) {                    // p=[kp,+inf)
        return F_14(p,k0-p);
    },  kp    )(x);
    res +=
    remap([&](double p, double pd) {                    // p=(-inf,km]
        return F_14(p,k0-p);
    },  km    )(x);
  }

  if ( isinf(res)||isnan(res) ) { return 0.;}
  else { return res/(k); }
}
