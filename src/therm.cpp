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
  res = gsl_sf_zeta_int(n+2)/( 2.*SQR(M_PI) );
  for (int i=n+1;i>0;i--) res*=(double)i;
  res*= (double)s_ + (s_<0?pow(2.,-n-1):0.);
  return res;
}
