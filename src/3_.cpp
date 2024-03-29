
#include "core.hh"
#include "quad.hh"
#include "map.hh"

struct H0 {
  // integrand for function H^[0]_m
  int nu, sA, sB;
  double muA = 0., muB = 0.;
  double operator ()(double p, double);
  double operator ()(double p);
  double val(int _nu, int _sA, int _sB);
}; 

struct H1 {
  // integrand for function H^[1]_m
  // O(eps) terms
  int nu, sA, sB;
  double muA = 0., muB = 0.;
  double operator ()(double p, double del);
  double operator ()(double p);
  double val(int _nu, int _sA, int _sB);
};

struct G0 {
  // integrand for G^[0]_m
  int nu, sA, sB;
  double muA = 0., muB = 0.;
  double operator ()(double p, double del);
  double operator ()(double p);
  double val(int _nu, int _sA, int _sB);
};

/*--------------------------------------------------------------------*/

struct rho11011 : Master {
  double integrand(double,double) { return 0.; };
  double eval();

  H0 h0, hh0;
  H1 h1, hh1;
  G0 g0, gg0;

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

  h0.muA = MOT1; h0.muB = MOT4;
  h1.muA = MOT1; h1.muB = MOT4;
  g0.muA = MOT1; g0.muB = MOT4;

  hh0.muA = MOT2; hh0.muB = MOT5;
  hh1.muA = MOT2; hh1.muB = MOT5;
  gg0.muA = MOT2; gg0.muB = MOT5;

  double a1=I(0,(this->s)[1],MOT1), b1=I(2,(this->s)[1],MOT1), // tadpole ints
         a2=I(0,(this->s)[2],MOT2), b2=I(2,(this->s)[2],MOT2),
         a4=I(0,(this->s)[4],MOT4), b4=I(2,(this->s)[4],MOT4),
         a5=I(0,(this->s)[5],MOT5), b5=I(2,(this->s)[5],MOT5);

  if ( m==0 && n==0 ) { // (0)
    (this->OPE).T0 = -2.;
    (this->OPE).T2 = ( a1+a2+a4+a5 )*.25*OOFP/K2;
    (this->OPE).T4 = ( b1+b2+b4+b5 )*(k0*k0+k*k/3.)/CUBE(K2)*OOFP;
  } else
  if ( m==1 && n==0 ) { // (1)
    (this->OPE).T0 = -1.;
    (this->OPE).T2 = ( a2+2.*a4+a5 )*.125*OOFP*k0/K2;
    (this->OPE).T4 = ( +(b1-b4)*K2 
                       +(b2+2.*b4+b5)*(k0*k0+k*k/3.) )*k0/CUBE(K2)*.5*OOFP;
  } else
  if ( m==0 && n==1 ) { // (0,1)
    (this->OPE).T0 = -1.;
    (this->OPE).T2 = ( a1+2.*a5+a4 )*.125*OOFP*k0/K2;
    (this->OPE).T4 = ( +(b2-b5)*K2 
                       +(b1+2.*b5+b4)*(k0*k0+k*k/3.) )*k0/CUBE(K2)*.5*OOFP;
  } else
  if ( m==1 && n==1 ) { // (1,1)
    (this->OPE).T0 = -.5;
    (this->OPE).T2 = ( a4+a5 )*.125*OOFP*SQR(k0)/K2;
    (this->OPE).T4 = ( -(b4+b5)*K2 
                       +(b4+b5)*(k0*k0+k*k/3.) )*SQR(k0)/CUBE(K2)*.5*OOFP;
  } else {
    cerr << "Case: (m,n)=("<< m << ","<<n<<") out of bounds!\n";
    return 0.;
  }
  //return (this->OPE)();

  h0m = h0.val(m,s[1],s[4]);  h0n = hh0.val(n,s[2],s[5]);
  h1m = h1.val(m,s[1],s[4]);  h1n = hh1.val(n,s[2],s[5]);
  g0m = g0.val(m,s[1],s[4]);  g0n = gg0.val(n,s[2],s[5]);
  g1m = pow(.5,m);            g1n = pow(.5,n);

  res += -(2.*g1m-g0m)*h0n+g1m*h1n;
  res += -(2.*g1n-g0n)*h0m+g1n*h1m;

  return -res*.25*CUBE(OOFP)*pow(k0,m+n);
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
        return pow(p/k0,nu)*f(p-muA,sA)*f(k0-p-muB,sB);
    }, km,kp )(x);
  } else {
    res =
    remap([&](double p, double del) {                 // p=[kp,+inf)
        return pow(p/k0,nu)*f(p-muA,sA)*f(k0-p-muB,sB);
    }, kp )(x);
    res +=
    remap([&](double p, double del) {                 // p=(-inf,km]
        return pow(p/k0,nu)*f(p-muA,sA)*f(k0-p-muB,sB);
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

  double epsabs = 1e-7, epsrel = 1e-8;
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
        temp  = pow(p/k0,nu)*f(p-muA,sA)*f(k0-p-muB,sB)  ;
        pm    = (fabs(p-km)<.1*k) ? del : (p-km) ;    // p=km (!)
        pp    = (fabs(p-kp)<.1*k) ? del : (p-kp) ;    //  =kp (!)
        temp *= lga(k*k/(pp*pm)); return temp;
    }, km,kp )(x);
  } else {
    res =
    remap([&](double p, double del) {                 // p=[kp,+inf)
        double temp;
        temp  = pow(p/k0,nu)*f(p-muA,sA)*f(k0-p-muB,sB)  ;
        pm    = (p-km) ;
        pp    = (fabs(p-kp)<.1*k) ? del : (p-kp) ;    // p=kp (!)
        temp *= lga(k*k/(pp*pm)); return temp;
    }, kp )(x);
    res +=
    remap([&](double p, double del) {                 // p=(-inf,km]
        double temp;
        temp  = pow(p/k0,nu)*f(p-muA,sA)*f(k0-p-muB,sB)  ;
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

  double epsabs = 1e-4, epsrel = 1e-5;
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
      //lg1 = lga(pm/pp);
      //lg2 = lga(SQR(kp/km))+lg1;
      lg1 = lga(pm*kp/(pp*km));

      pp  = p+kp; pm = p+km;
      if (km<0) pm  = fabs(pm)<1e-1*(-km) ? del : pm; // p=-km (!)
      //lg2-= lga(pp/pm);
      //lg1+= lga(pp/pm);
      lg2 = lga(pp*km/(pm*kp));

      if (nu==0) temp = ( ( f(p-muA,sA)+f(p-muB,sB) )*lg1
                        + ( f(p+muA,sA)+f(p+muB,sB) )*lg2 );
      if (nu==1) temp = ( f(p-muB,sB)*lg1 + (p/k0)*( f(p-muA,sA)-f(p-muB,sB) )*lg1
                        + f(p+muB,sB)*lg2 - (p/k0)*( f(p+muA,sA)-f(p+muB,sB) )*lg2 );
      if (nu==2) temp = ( SQR(1.-p/k0)*f(p-muB,sB)*lg1 + SQR(+p/k0)*f(p-muA,sA)*lg1
                        + SQR(1.+p/k0)*f(p+muB,sB)*lg2 + SQR(-p/k0)*f(p+muA,sA)*lg2 );
      return temp;
  }, 0,fabs(km) )(x);                                 //*/

  _B =
  remap([&](double p, double del) {                   // p=[|km|,kp]
      double temp=0., lg1, lg2;

      pp  = p-kp; pm = p-km; 
      if (km>0) pm = fabs(pm)<1e-1*k ? del : pm;      // p=+km (!)
      pp = fabs(pp)<1e-1*(kp-fabs(km)) ? del : pp;    // p= kp (!)
      //lg1 = lga(pm/pp);
      //lg2 = lga(SQR(kp/km))+lg1;
      lg1 = lga(pm*kp/(pp*km));

      pp  = p+kp; pm = p+km;
      if (km<0) pm = fabs(pm)<1e-1*k0 ? del : pm;     // p=-km (!)
      //lg2-= lga(pp/pm);
      //lg1+= lga(pp/pm);
      lg2 = lga(pp*km/(pm*kp));

      //if (nu==0) temp = ( f(p-muA,sA)+f(p-muB,sB) )*lg1;
      if (nu==0) temp = ( ( f(p-muA,sA)+f(p-muB,sB) )*lg1
                        + ( f(p+muA,sA)+f(p+muB,sB) )*lg2 );
      //if (nu==1) temp = f(p-muB,sB)*lg1 + (p/k0)*( f(p-muA,sA)-f(p-muB,sB) )*lg2;
      if (nu==1) temp = ( f(p-muB,sB)*lg1 + (p/k0)*( f(p-muA,sA)-f(p-muB,sB) )*lg1
                        + f(p+muB,sB)*lg2 - (p/k0)*( f(p+muA,sA)-f(p+muB,sB) )*lg2 );
      if (nu==2) temp = ( SQR(1.-p/k0)*f(p-muB,sB)*lg1 + SQR(+p/k0)*f(p-muA,sA)*lg1
                        + SQR(1.+p/k0)*f(p+muB,sB)*lg2 + SQR(-p/k0)*f(p+muA,sA)*lg2 );
      return temp;
  }, fabs(km),kp )(x);                                //*/

  _C =
  remap([&](double p, double del) {                   // p=[kp,+inf)
      double temp=0., lg1, lg2;

      pp  = p-kp; pm = p-km;
      pp = fabs(pp)<1e-1*k0 ? del : pp;               // p=+kp (!)
      //lg1 = lga(pm/pp);
      //lg2 = lga(SQR(kp/km))+lg1;
      lg1 = lga(pm*kp/(pp*km));

      pp  = p+kp; pm = p+km;
      //lg2-= lga(pp/pm);
      //lg1+= lga(pp/pm);
      lg2 = lga(pp*km/(pm*kp));

      //if (nu==0) temp = ( f(p-muA,sA)+f(p-muB,sB) )*lg1;
      if (nu==0) temp = ( ( f(p-muA,sA)+f(p-muB,sB) )*lg1
                        + ( f(p+muA,sA)+f(p+muB,sB) )*lg2 );
      //if (nu==1) temp = f(p-muA,sB)*lg1 + (p/k0)*( f(p-muA,sA)-f(p-muB,sB) )*lg2;
      if (nu==1) temp = ( f(p-muB,sB)*lg1 + (p/k0)*( f(p-muA,sA)-f(p-muB,sB) )*lg1
                        + f(p+muB,sB)*lg2 - (p/k0)*( f(p+muA,sA)-f(p+muB,sB) )*lg2 );
      if (nu==2) temp = ( SQR(1.-p/k0)*f(p-muB,sB)*lg1 + SQR(+p/k0)*f(p-muA,sA)*lg1
                        + SQR(1.+p/k0)*f(p+muB,sB)*lg2 + SQR(-p/k0)*f(p+muA,sA)*lg2 );
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
  //cout << nu << endl;

  double epsabs = 1e-5, epsrel = 1e-4;
  size_t limit = 1e5;

  quad wsp(limit);

  auto outer = make_gsl_function( [&](double x) { return (*this)(x); } );
  gsl_integration_qag(  outer, .0,1., epsabs, epsrel, 
                        limit, 6, wsp, &res, &err  );

  return res;
};

/*--------------------------------------------------------------------*/

double chi(int nu, int sA, int sB) {
  //cout << "((TEST))" << endl;
  G0 g0;
  g0.muA = 0.; g0.muB = 0.;
  double res = g0.val(nu,sA,sB);
  //double res = 1.;
  return res/SQR(4.*M_PI);
}

