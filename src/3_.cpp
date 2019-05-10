
#include "core.hh"
#include "quad.hh"
#include "map.hh"

struct H0 {
  // integrand for function H^[0]_m
  int nu, sA, sB;
  double operator ()(double p, double);
  double operator ()(double p);
  double val(int _nu, int _sA, int _sB);
}; 

struct H1 {
  // integrand for function H^[1]_m
  // O(eps) terms
  int nu, sA, sB;
  double operator ()(double p, double del);
  double operator ()(double p);
  double val(int _nu, int _sA, int _sB);
};

struct G0 {
  // integrand for G^[0]_m
  int nu, sA, sB;
  double operator ()(double p, double del);
  double operator ()(double p);
  double val(int _nu, int _sA, int _sB);
};

/*--------------------------------------------------------------------*/

struct rho11011 : Master {

  double eval();

  H0 h0;
  H1 h1;
  G0 g0;
  double h0m, h0n, h1m, h1n, g0m, g0n, g1m, g1n;
  map *S; 
  rho11011(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=3; }
};
// function for MAIN
Master* _11011(int m, int n, int s[3]) {
  Master *R =  new rho11011(m,n,s); return R;
}

double rho11011::eval() 
{
  double res=0.;

  double a1=I(0,(this->s)[1]), b1=I(2,(this->s)[1]), // tadpole ints
         a2=I(0,(this->s)[2]), b2=I(2,(this->s)[2]),
         a4=I(0,(this->s)[4]), b4=I(2,(this->s)[4]),
         a5=I(0,(this->s)[5]), b5=I(2,(this->s)[5]);

  if ( m==0 && n==0 ) { // (0)
    (this->OPE).T0 = -2.*K2;
    (this->OPE).T2 = ( a1+a2+a4+a5 )*.25*OOFP/K2;
    (this->OPE).T4 = ( b1+b2+b4+b5 )*(k0*k0+k*k/3.)/CUBE(K2)*OOFP;
  } else
  if ( m==1 && n==0 ) { // (1)
    (this->OPE).T0 = -K2;
    (this->OPE).T2 = ( a2+2.*a4+a5 )*.125*OOFP/K2;
    (this->OPE).T4 = ( -(b1-b4)*K2 
                       +(b2+2.*b4+b5)*(k0*k0+k*k/3.) )/CUBE(K2)*.5*OOFP;
  } else
  if ( m==0 && n==1 ) { // (0,1)
    (this->OPE).T0 = -K2;
    (this->OPE).T2 = ( a1+2.*a5+a4 )*.125*OOFP/K2;
    (this->OPE).T4 = ( +(b2-b5)*K2 
                       +(b1+2.*b5+b4)*(k0*k0+k*k/3.) )/CUBE(K2)*.5*OOFP;
  } else
  if ( m==1 && n==1 ) { // (1,1)
    (this->OPE).T0 = -.5*K2;
    (this->OPE).T2 = ( a4+a5 )*.125*OOFP/K2;
    (this->OPE).T4 = ( -(b4+b5)*K2 
                       +(b4+b5)*(k0*k0+k*k/3.) )/CUBE(K2)*.5*OOFP;
  } else {
    cerr << "Case: (m,n)=("<< m << ","<<n<<") out of bounds!\n";
    return 0.;
  }

  h0m = h0.val(m,s[1],s[4]);  h0n = h0.val(n,s[2],s[5]);
  h1m = h1.val(m,s[1],s[4]);  h1n = h1.val(n,s[2],s[5]);
  g0m = g0.val(m,s[1],s[4]);  g0n = g0.val(n,s[2],s[5]);
  g1m = pow(.5,m);            g1n = pow(.5,n);

  res += -(2.*g1m-g0m)*h0n+g1m*h1n;
  res += -(2.*g1n-g0n)*h0m+g1n*h1m;

  return -res*.25*CUBE(OOFP);
}

/*--------------------------------------------------------------------*/
// here, the helper functions:


// H0
double H0::operator ()(double x)
{
  double res;

  if (k0>k) {
    res =
    remap([&](double p, double del) {                 // p=[km,kp]
        return pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
    }, km,kp )(x);
  } else {
    res =
    remap([&](double p, double del) {                 // p=[kp,+inf)
        return pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
    }, kp )(x);
    res +=
    remap([&](double p, double del) {                 // p=(-inf,km]
        return pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
    }, km )(x);
  }
  return (((double)sA*sB)*exp(k0)-1.)*res/k*sgn(km);  //*/
}

double H0::val(int _nu, int _sA, int _sB)
{
  nu=_nu;
  sA=_sA;
  sB=_sB;
  double res=0.,err;

  double epsabs = 1e-3, epsrel = 1e-6;
  size_t limit = 1e5;

  quad wsp(limit);
  auto outer = make_gsl_function( [&](double x) {
      return (*this)(x);
      } );
  gsl_integration_qag(  outer, .0,1., epsabs, epsrel, 
                        limit, 6, wsp, &res, &err  );
  return res;
}

// H1
double H1::operator ()(double x)
{
  double res, pp, pm; // pi_+=(p-k_+), pi_-=...

  if (k0>k) {
    res =
    remap([&](double p, double del) {                 // p=[km,kp]
        double temp;
        temp  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB)  ;
        pm    = (fabs(p-km)<.1*k) ? del : (p-km) ;    // p=km (!)
        pp    = (fabs(p-kp)<.1*k) ? del : (p-kp) ;    //  =kp (!)
        temp *= lga(k*k/(pp*pm)); return temp;
    }, km,kp )(x);
  } else {
    res =
    remap([&](double p, double del) {                 // p=[kp,+inf)
        double temp;
        temp  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB)  ;
        pm    = (p-km) ;
        pp    = (fabs(p-kp)<.1*k) ? del : (p-kp) ;    // p=kp (!)
        temp *= lga(k*k/(pp*pm)); return temp;
    }, kp )(x);
    res +=
    remap([&](double p, double del) {                 // p=(-inf,km]
        double temp;
        temp  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB)  ;
        pm    = (fabs(p-km)<.1*k) ? del : (p-km) ;    // p=km (!)
        pp    = (p-kp) ;
        temp *= lga(k*k/(pp*pm)); return temp;
    }, km )(x);
  }

  return (((double)sA*sB)*exp(k0)-1.)*res/k*sgn(km);  //*/
}

double H1::val(int _nu, int _sA, int _sB)
{
  nu=_nu;
  sA=_sA;
  sB=_sB;
  double res, err;

  double epsabs = 1e-2, epsrel = 1e-2;
  size_t limit = 1e5;

  quad wsp(limit);

  auto outer = make_gsl_function( [&](double x) { return (*this)(x); } );
  gsl_integration_qag(  outer, .0,1., epsabs, epsrel,
                        limit, 6, wsp, &res, &err  ); // 61 pt. Gauss-Konrod
  return res;
}

// G0 (nu=0,1)
double G0::operator ()(double x) 
{
  // I construct p(x) by subdividing the interval [0,+inf) into
  //
  //   [0,|km|] U [|km|,kp] U [kp,+inf)
  //   A          B           C
  //

  double pp, pm;                                      // = p \pm k_\pm , resp.

  double _A, _B, _C;

  _A =
  remap([&](double p, double del) {                   // p=[0,|km|]
      double temp=0., lg1, lg2;

      pp  = p-kp; pm = p-km;
      if (km>0) pm  = fabs(pm)<1e-1*(+km) ? del : pm; // p=+km (!)
      lg1 = lga(pm/pp);
      lg2 = lga(SQR(kp/km))+lg1;

      pp  = p+kp; pm = p+km;
      if (km<0) pm  = fabs(pm)<1e-1*(-km) ? del : pm; // p=-km (!)
      lg2-= lga(pp/pm);
      lg1+= lga(pp/pm);

      if (nu==0) temp = ( f(p,sA)+f(p,sB) )*lg1;
      if (nu==1) temp = f(p,sB)*lg1 + (p/k0)*( f(p,sA)-f(p,sB) )*lg2;
      return temp;
  }, 0,fabs(km) )(x);                                 //*/

  _B =
  remap([&](double p, double del) {                   // p=[|km|,kp]
      double temp=0., lg1, lg2;

      pp  = p-kp; pm = p-km; 
      if (km>0) pm = fabs(pm)<1e-1*k ? del : pm;      // p=+km (!)
      pp = fabs(pp)<1e-1*(kp-fabs(km)) ? del : pp;    // p= kp (!)
      lg1 = lga(pm/pp);
      lg2 = lga(SQR(kp/km))+lg1;

      pp  = p+kp; pm = p+km;
      if (km<0) pm = fabs(pm)<1e-1*k0 ? del : pm;     // p=-km (!)
      lg2-= lga(pp/pm);
      lg1+= lga(pp/pm);

      if (nu==0) temp = ( f(p,sA)+f(p,sB) )*lg1;
      if (nu==1) temp = f(p,sB)*lg1 + (p/k0)*( f(p,sA)-f(p,sB) )*lg2;
      return temp;
  }, fabs(km),kp )(x);                                //*/

  _C =
  remap([&](double p, double del) {                   // p=[kp,+inf)
      double temp=0., lg1, lg2;

      pp  = p-kp; pm = p-km;
      pp = fabs(pp)<1e-1*k0 ? del : pp;               // p=+kp (!)
      lg1 = lga(pm/pp);
      lg2 = lga(SQR(kp/km))+lg1;

      pp  = p+kp; pm = p+km;
      lg2-= lga(pp/pm);
      lg1+= lga(pp/pm);

      if (nu==0) temp = ( f(p,sA)+f(p,sB) )*lg1;
      if (nu==1) temp = f(p,sB)*lg1 + (p/k0)*( f(p,sA)-f(p,sB) )*lg2;
      return temp;
  }, kp )(x);                                         //*/

  return ( _A + _B + _C )/k;
};

double G0::val(int _nu, int _sA, int _sB)
{
  nu=_nu;
  sA=_sA;
  sB=_sB;
  double res, err;

  double epsabs = 1e-2, epsrel = 1e-2;
  size_t limit = 1e5;

  quad wsp(limit);

  auto outer = make_gsl_function( [&](double x) { return (*this)(x); } );
  gsl_integration_qag(  outer, .0,1., epsabs, epsrel, 
                        limit, 6, wsp, &res, &err  );

  return res;
};
