
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
  double h0m, h0n, h1m, h1n, g0m, g0n, g1m, g1n;
  map *S; 
  rho11011(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=3; }
};
// function for MAIN
Master* _11011(int m, int n, int s[3]) {
  Master *R =  new rho11011(m,n,s); return R;
}

double rho11011::eval() {
  double res=0.;
  h0m=h0.val(m,s[1],s[4]); h0n=h0.val(n,s[2],s[5]);
  h1m=h1.val(m,s[1],s[4]); h1n=h1.val(n,s[2],s[5]);
  g0m=g0.val(m,s[1],s[4]); g0n=g0.val(n,s[2],s[5]);
  g1m=pow(.5,m);        g1n=pow(.5,n);
  res += -(2.*g1m-g0m)*h0n+g1m*h1n;
  res += -(2.*g1n-g0n)*h0m+g1n*h1m;
  return -res*.25*pow(OOFP,3);
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
  double res, pp, pm; // pi_+=(p-k_+), pi_-=...
  res  = pow(p/k0,nu)*f(p,sA)*f(k0-p,sB)  ;
  pp   = (fabs(p-km)<.1*k) ? del : (p-km) ;
  pm   = (fabs(p-kp)<.1*k) ? del : (p-kp) ;
  res *= lga(k*k/(pp*pm));
  //if      (fabs(p-km)<.1*k) { res *= lga(k*k/( (p-kp)*del )); }
  //else if (fabs(p-kp)<.1*k) { res *= lga(k*k/( (p-km)*del )); }
  //else                      { return (*this)(p); }
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
  if (nu==0) {
    res  = f(p,sA)+f(p,sB);
    res *= lga( (p-km)*(p+kp)/((p-kp)*(km+p)) ); 
  } else 
  if (nu==1) {
    res  = f(p,sB)*lga( (p-km)*(p+kp)/((p-kp)*(km+p)) ); 
    res += (p/k0)*( f(p,sA)-f(p,sB) )*lga( kp*kp*(km*km-p*p) /( km*km*(kp*kp-p*p) ) );
  }
  return res/k;
};

double G0::operator ()(double p, double del) 
  // for singularities, del parametrises p-s
  // TODO: tify up if-elses
{
  double res, h;
  h = .1*fabs(kp-fabs(km));
  // singularities @ (s-h,s+h) are assisted:
  // ( s={abs(km),kp} for p>0 )
  if (nu==0) {
    res  = f(p,sA)+f(p,sB);
    if (fabs(p-kp)<h) // close to p=kp singularity
      res *= lga( (p-km)*(p+kp)/(del*(km+p)) );
    else if (( (fabs(p-fabs(km)))<h&&(p>fabs(km)/2.) )) {
      if (km>0.)  // close to p=km singularity
        res *= lga( del*(p+kp)/((p-kp)*(km+p)) );
      else // close to p=-km
        res *= lga( (p-km)*(p+kp)/((p-kp)*del) ); 
    }
    else // "safe" (or p=0 singularity)
      return (*this)(p); 
    res /= k; return res;
  } else
  if (nu==1) {
    res  = 0.;
    if (fabs(p-kp)<h) {// close to p=kp singularity, del=(kp-p)
      res += f(p,sB)*lga( (p-km)*(p+kp)/(del*(km+p)) );
      res += (p/k0)*( f(p,sA)-f(p,sB) )*lga( kp*kp*(p-km)*(p+km)/(km*km*del*(kp+p)) );
    }
    else if (( (fabs(p-fabs(km)))<h&&(p>fabs(km)/2.) )) {
      if (km>0.) { // close to p=km singularity, del=(km-p)
        res += f(p,sB)*lga( del*(p+kp)/((p-kp)*(km+p)) );
        res += (p/k0)*( f(p,sA)-f(p,sB) )*lga( kp*kp*del*(p+km)/(km*km*(p-kp)*(kp+p)) );
      }
      else {// close to p=-km, del=(km+p)
        res += f(p,sB)*lga( (p-km)*(p+kp)/((p-kp)*del) );
        res += (p/k0)*( f(p,sA)-f(p,sB) )*lga( kp*kp*(p-km)*del/(km*km*(p-kp)*(p+kp)) ); 
      }
    }
    else // "safe" (or p=0 singularity)
      return (*this)(p); 
    res /= k; return res;
  } 
  return 0;
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
