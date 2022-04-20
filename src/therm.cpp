#include "core.hh"
#include <math.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>

/*
 *  distributions functions
 *
 *--------------------------------------------------------------------*/

double nB( double x ) {
  // Bose-Einstein
  if (x>20.) {double e=exp(-x); return e/(1.-e);}
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

double I(int n, int s_, double mu_ = 0.) {
  // tadpole, needed for OPEs
  double res;
  res = gsl_sf_zeta_int(n+2)/( 2.*SQR(M_PI) );
  for (int i=n+1;i>0;i--) res*=(double)i;
  res*= (double)s_ + (s_<0?pow(2.,-n-1):0.);
  if (s_<0) { 
    if (n==0) res *= 1.+3.*SQR(4.*OOFP*mu_);
    if (n==2) res *= 1.+60.*SQR(4.*OOFP*mu_)/7.+30.*SQR(SQR(4.*OOFP*mu_))/7.;
  }
  return res;
}

/*
 *  useful moments
 *
 *--------------------------------------------------------------------*/
#define Ep exp(-fabs(kp))
#define Em exp(-fabs(km))

double Li2(double z) { return gsl_sf_dilog(z); }

double Li3(double z) { // z=(-\infty,1]
                       // NB: z>1 is not included!!!
  double res = 0., temp, Z3 = gsl_sf_zeta(3), lz = log(fabs(z));

  if (z<-1.)        { return Li3(1./z) - lz*( SQR(lz) + SQR(M_PI) )/6. ; }
  else if (z==1.)   { return      Z3; }
  else if (z==-1.)  { return -.75*Z3; }
  else if (z<0.)    { return .25*Li3(SQR(z))-Li3(-z); }
  else if (z<.25)   {
    temp = z;
    for (int i=1;i<10;i++) {
      res  += temp/CUBE((double)i);
      temp *= z;
    }
  }
  else if (z>.25)   {
    res   = Z3 + SQR(M_PI)*lz/6.+.5*SQR(lz)*(1.5-log(-lz));
    temp  = 2.;
    for (int i=3;i<10;i++) {
      temp *= (double)i;
      res+= gsl_sf_zeta(3-i)*pow(lz,i)/temp;
    }
  }
  return res;
}

double l1(int s, double x) {
  return lga(1.-((double)s)*exp(-x));
}

double l2(int s, double x) {
  return Li2(((double)s)*exp(-x));
}

double l3(int s, double x) {
  return Li3(((double)s)*exp(-x));
}


double psi0(int sA, int sB, double muA = 0.) {
  double res  = (km>0) ? 1. : 0.;
  //double exA = exp(-muA), exB = exp(+muA), exC = exp(-sgn(km)*muA);
  //res += lga( (1.-((double)sA)*Ep*exA)/(1.-((double)sA)*Em*exC) )/k;
  //res += lga( (1.-((double)sB)*Ep*exB)/(1.-((double)sB)*Em/exC) )/k;
  double sigma = sgn(km);
  res += ( l1(sA,kp-muA) - l1(sA,fabs(km)-sigma*muA) )/k;
  res += ( l1(sB,kp+muA) - l1(sB,fabs(km)+sigma*muA) )/k;
  return res;
}

double psi1(int sA, int sB, double muA = 0.) {
  double res  = (km>0) ? .5 : 0.;
  //double ap = ((double)sA)*Ep,
  //       am = ((double)sA)*Em,
  //       bp = ((double)sB)*Ep,
  //       bm = ((double)sB)*Em;
  //res += k0*lga( (1.-bp)/(1.-bm) )/(k*k0);
  //res += kp*lga( (1.-ap)/(1.-bp) )/(k*k0);
  //res -= km*lga( (1.-am)/(1.-bm) )/(k*k0);
  //res -= ( Li2(ap)+Li2(am)*sgn(km)-Li2(bp)-Li2(bm)*sgn(km) )/(k*k0);
  double sigma = sgn(km);
  res += ( kp*l1(sA,kp-muA) - km*l1(sA,fabs(km)-sigma*muA) )/(k*k0);
  res += ( km*l1(sB,kp+muA) - kp*l1(sA,fabs(km)+sigma*muA) )/(k*k0);
  res += ( -l2(sA,kp-muA) +sigma*l2(sA,fabs(km)-sigma*muA) )/(k*k0);
  res += ( +l2(sB,kp+muA) -sigma*l2(sB,fabs(km)+sigma*muA) )/(k*k0);
  return res;
}

double psi2(int sA, int sB, double muA = 0.) {
  double res  = (km>0) ? .25 : 0.;
  //double ap = ((double)sA)*Ep,
  //       am = ((double)sA)*Em,
  //       bp = ((double)sB)*Ep,
  //       bm = ((double)sB)*Em;
  //res *= ( 1.+SQR(k/k0)/3. );
  //res += SQR(kp/k0)*lga( (1.-ap)/(1.-bm) )/k;
  //res -= SQR(km/k0)*lga( (1.-am)/(1.-bp) )/k;
  //res -= 2.*kp*( Li2(ap)+Li2(bm)*sgn(km) )/(k*SQR(k0));
  //res += 2.*km*( Li2(bp)+Li2(am)*sgn(km) )/(k*SQR(k0));
  //res -= 2.*( Li3(ap)-Li3(am)+Li3(bp)-Li3(bm) )/(k*SQR(k0));
  double sigma = sgn(km);
  res *= ( 1.+SQR(k/k0)/3. );
  res += ( SQR(kp)*l1(sA,kp-muA) - SQR(km)*l1(sA,fabs(km)-sigma*muA) )/(k*SQR(k0));
  res += ( SQR(km)*l1(sB,kp+muA) - SQR(kp)*l1(sB,fabs(km)+sigma*muA) )/(k*SQR(k0));
  res += 2.*( - kp*l2(sA,kp-muA) + sigma*km*l2(sA,fabs(km)-sigma*muA) )/(k*SQR(k0));
  res += 2.*( + km*l2(sB,kp+muA) - sigma*kp*l2(sB,fabs(km)+sigma*muA) )/(k*SQR(k0));
  res += 2.*( - l3(sA,kp-muA) + l3(sA,fabs(km)-sigma*muA) )/(k*SQR(k0));
  res += 2.*( - l3(sB,kp+muA) + l3(sB,fabs(km)+sigma*muA) )/(k*SQR(k0));
  return res;
}
