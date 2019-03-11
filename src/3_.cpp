
#include "core.hh"
#include "trapezoid.hh"
#include "map.hh"
/*
template<typename F>
double in(F *func) {
  map *S; 
  double res=0.;
    if (k0>k) { 
      S = new Finite<F>(*func,km,kp);
      integrate<map> i(*S);
      res += go(i);
    } else {
      S = new SemiInf<F>(*func,km);
      integrate<map> i1(*S);
      res += go(i1);
      S = new SemiInf<F>(*func,kp);
      integrate<map> i2(*S);
      res += go(i2);
    }
  return res;
}*/

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
  map *S; 
  rho11011(int _m, int _n, int _s[3]) : Master(_m,_n,_s) {
   // g0m = new G0(m,s[1],s[4]);
  }
};
// function for MAIN
Master* _11011(int m, int n, int s[3]) {
  Master *R =  new rho11011(m,n,s); return R;
}

double rho11011::eval() {
  double res=0.;
  res += h0.val(0,s[1],s[4]);
  return res;
}

/*--------------------------------------------------------------------*/
//
// here the helper functions:

double H0::operator ()(double p)
{
  double res;
  res  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
  return (((double)sA*sB)*exp(k0)-1.)*res/k;
}

double H0::operator ()(double p, double)
{
  return (*this)(p);
}

double H0::val(int _nu, int _sA, int _sB) {
  nu=_nu;
  sA=_sA;
  sB=_sB;
  double res=0.;
  if (k0>k) { 
    Finite<H0> S(*this,km,kp);
    integrate<map> i(S);
    res += go(i);
  } else {
    SemiInf<H0> S(*this,km);
    integrate<map> i1(S);
    SemiInf<H0> SS(*this,kp);
    integrate<map> i2(SS);
    res += go(i1)-go(i2);
  }
  return res;
}


double H1::operator ()(double p)
{
  double res;
  res  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
  res *= lga(k*k/( (p-kp)*(p-km) ));
  return (((double)sA*sB)*exp(k0)-1.)*res/k;
}

double H1::operator ()(double p, double del)
{
  double res;
  res  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB);
  if      (fabs(p-km)<.1*k) { res *= lga(k*k/( (p-kp)*del )); }
  else if (fabs(p-kp)<.1*k) { res *= lga(k*k/( (p-km)*del )); }
  else                      { return (*this)(p); }
  return (((double)sA*sB)*exp(k0)-1.)*res/k;
}

double H1::val(int _nu, int _sA, int _sB) {
  nu=_nu;
  sA=_sA;
  sB=_sB;
  double res=0.;
  if (k0>k) { 
    Finite<H1> S(*this,km,kp);
    integrate<map> i(S);
    res += go(i);
  } else {
    SemiInf<H1> S(*this,km);
    integrate<map> i1(S);
    SemiInf<H1> SS(*this,kp);
    integrate<map> i2(SS);
    res += go(i1)-go(i2);
  }
  return res;
}

double G0::operator ()(double p)
  // raw integrand, function of p
{
  double res;
  res  = f(p,sA)+f(p,sB);
  res *= lga( (p-km)*(p+kp)/((p-kp)*(km+p)) ); 
  return res/k;
};

double G0::operator ()(double p, double del) 
  // for singularities, del parametrises p-s
{
  double res, h;
  res  = f(p,sA)+f(p,sB);
  res /= k;
  h = .1*fabs(kp-fabs(km));
  // singularities @ (s-h,s+h) are assisted:
  // ( s={abs(km),kp} for p>0 )
  if (fabs(p-kp)<h) // close to p=kp singularity
    return res*lga( (p-km)*(p+kp)/(del*(km+p)) );
  else if (( (fabs(p-fabs(km)))<h&&(p>fabs(km)/2.) )) {
    if (km>0.)
      // close to p=km singularity
      return res*lga( del*(p+kp)/((p-kp)*(km+p)) );
    else // close to p=-km
      return res*lga( (p-km)*(p+kp)/((p-kp)*del) ); }
  else // "safe" (or p=0 singularity)
    return (*this)(p); 
};

double G0::val(int _nu, int _sA, int _sB)
  // evaluate the integral
{
  nu=_nu;
  sA=_sA;
  sB=_sB;
  //if (k0>k) { // split ip the integration range
    //: [0,km]+[km,kp]+[kp,+inf)
    Finite<G0> S(*this,0.,fabs(km));
    integrate<map> i1(S);
    Finite<G0> SS(*this,fabs(km),kp);
    integrate<map> i2(SS);
    SemiInf<G0> SSS(*this,kp);
    integrate<map> i3(SSS);
    return go(i1)+go(i2)+go(i3);
};
