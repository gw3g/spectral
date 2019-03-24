#include "core.hh"
#include <math.h>
#include <gsl/gsl_sf_zeta.h>

/*
 *  distributions functions
 *
 *--------------------------------------------------------------------*/

double nB( double x ) {
  // Bose-Einstein
  if (x>10.) {double e=exp(-x); return e/(1.-e);}
  else return 1./expm1(x) ;
};

double nF( double x ) {
  // Fermi-Dirac
  if (x>10) {double e=exp(-x); return e/(1.+e);}
  else return 1./(exp(x)+1.);
};

double f(double e, int X) {
  switch(X) {       case (-1): return -nF(e);
                    case (+1): return +nB(e); }
  return exp(-e);
};

double bf(double e, int X) {
  switch(X) {       case (-1): return 1.-nF(e);
                    case (+1): return 1.+nB(e); }
  return exp(-e);
};

double I(int n, int s_) {
  // tadpole
  double res;
  res = gsl_sf_zeta_int(n+2)/(2.*M_PI*M_PI);
  for (int i=n+1;i>0;i--) res*=(double)i;
  res*= (double)s_ + (s_<0?pow(2.,-n-1):0.);
  return res;
}

double ff(double p, double l, int s1, int s4) {
  int s0=s1*s4;
  double fp=f(p,s1), fl=f(l,s4), f0=f(p+l,s0);
  return ((double)s0)*fp*fl/f0;
}

double fff(double p, double q, double r, int m, int n, int s[3]) {
  /*
   * [gain - loss] term, collected into
   * one factor via detailed balance.
   *
   */
  double res;
  double fp=f(p,s[0]), fq=f(q,s[1]), fr= f(r,s[2]);
  //res = 1. + fp + fq + fr + fp*fr + fq*fr + fp*fq;
  res = ( ((double) s[0]*s[1]*s[2])*exp(k0)-1. )*fp*fq*fr;
  //res = ( ((double) s0)*exp(k0)-1. )*f(p,sa)*f(q,sb)*f(r,s0*sa*sb);
  res*= pow(k0,-m-n)*pow(p,m)*pow(q,n) ;
  //if ( (p>0)&&(q>0)&&(r>0)&&(k0>k) ) {res-=1.;} // subtract vacuum
  return res;
}
