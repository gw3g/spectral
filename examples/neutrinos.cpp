#include "core.hh"
#include "quad.hh"
#include "map.hh"
#include "timer.hh"
#include <fstream>
#include <string>

double k0,k;
bool  CHEM_POT = false;
double MOT1, MOT2, MOT3, MOT4, MOT5;// mu over T
double tol = 1e-1; // absolute tolerance

using namespace std;

ofstream fout;
ifstream fin;
int elapsed; float percentage;

// for computing A,B,C,D from  2312.07015:

int Print_nu_rate(double,double,int,string,double);
//void QED_integrand_(double kk, double p_m, double p_p, double *_lo, double *_nlo);
void QED_integrand_(double kk, double p_m, double p_p, double *A, double *B, double *C, double *D);
double _A_integrand(double,double,double);
double _B_integrand(double,double,double);
double _C_integrand(double,double,double);
double _D_integrand(double,double,double);

double region_s(double,double (*)(double,double,double));
double region_t(double,double (*)(double,double,double));
double _D_t(double,double);
double _A_t(double);
double _B_t(double);

// real parts of 1-loop SPF:
double chi_V(double,double);
double chi_00(double,double);
double chi_T(double,double);
double chi_L(double,double);

int main(int argc, char *argv[]) {

  // Enter: k0, k
  if (argc>1) { k0=atof(argv[1]); k=atof(argv[2]); }

  /*double a1, a2;
  QED_integrand_(2.,.9,.1,&a1,&a2);
  cout << " pm = " << 4.3 << " , pp = " << 1. << endl;
  cout << " res = " << a1 << endl << endl;

  QED_integrand_(2.,-.9,.1,&a1,&a2);
  cout << " pm = " << -4.3 << " , pp = " << 1. << endl;
  cout << " res = " << a1 << endl;//*/

  //Print_nu_rate(.05,1.05,20,"QED_rho_1d",0.);
  Print_nu_rate(.5,15.,30,"QED_rho_rate",0.);
  //Print_nu_rate(.6,2.,1,"QED_rho_small_k_12",0.);
  //Print_nu_rate(.05,15.,5*60,"QED_rho_fix_Dt",0.);
  //
  //cout << " chi_V  = " <<  chi_V(-1.23,.179) << endl;
  //cout << " chi_00 = " << chi_00(-1.23,.179) << endl;

  /* for comparing with chi_test.dat from ML ... */
  /*
  double p0T, pT;

  cout << "# p/T" 
       << "  p0/T" 
       << "  chiV/T2" 
       << "  chi00/T2" 
       << "  chiT/T2" 
       << "  chiL/T2" 
       << endl;

  pT = .2; p0T = .1;
  cout << scientific
       << "  " <<  p0T
       << "  " <<  pT
       << "  " << chi_V(p0T,pT)
       << "  " << chi_00(p0T,pT)
       << "  " << chi_T(p0T,pT)
       << "  " << chi_L(p0T,pT)
       << endl;

  pT = 2.; p0T = .1;
  cout << scientific
       << "  " <<  p0T
       << "  " <<  pT
       << "  " << chi_V(p0T,pT)
       << "  " << chi_00(p0T,pT)
       << "  " << chi_T(p0T,pT)
       << "  " << chi_L(p0T,pT)
       << endl;

  pT = 2.; p0T = 1.;
  cout << scientific
       << "  " <<  p0T
       << "  " <<  pT
       << "  " << chi_V(p0T,pT)
       << "  " << chi_00(p0T,pT)
       << "  " << chi_T(p0T,pT)
       << "  " << chi_L(p0T,pT)
       << endl;

  pT = 2.; p0T = 3.;
  cout << scientific
       << "  " <<  p0T
       << "  " <<  pT
       << "  " << chi_V(p0T,pT)
       << "  " << chi_00(p0T,pT)
       << "  " << chi_T(p0T,pT)
       << "  " << chi_L(p0T,pT)
       << endl;//*/

  double aa,bb,cc,dd;
  //QED_integrand_(.5,-1.,.3,&aa,&bb,&cc,&dd);
  //QED_integrand_(.5,.7,1.3,&aa,&bb,&cc,&dd);
  //QED_integrand_(M_PI,-1.,.3,&aa,&bb,&cc,&dd);

  //cout << " A = " << aa << endl;
  //cout << " B = " << bb << endl;
  //cout << " C = " << cc << endl;
  //cout << " D = " << dd << endl;

  //cout << " A = " << _A_integrand(M_PI,-1.,.3) << endl;
  //cout << " B = " << _B_integrand(M_PI,-1.,.3) << endl;
  //cout << " C = " << _C_integrand(M_PI,-1.,.3) << endl;
  //cout << " D = " << _D_integrand(M_PI,-1.,.3) << endl;


/*
  double k_temp = 2.;
  cout << scientific << " A(s) = " << region_s(k_temp,_A_integrand) << endl;
  cout << scientific << " A(t) = " << region_t(k_temp,_A_integrand) << endl << endl;
  cout << scientific << " A(t)[imp] = " << _A_t(k_temp) << endl << endl;

  //cout << " C(s) = " << region_s(k_temp,_C_integrand) << endl;
  //cout << " C(t) = " << region_t(k_temp,_C_integrand) << endl << endl;

  cout << scientific << " B(s) = " << region_s(k_temp,_B_integrand) << endl;
  cout << scientific << " B(t) = " << region_t(k_temp,_B_integrand) << endl << endl;
  cout << scientific << " B(t)[imp] = " << _B_t(k_temp) << endl << endl;

  //cout << " D(s) = " << region_s(k_temp,_D_integrand) << endl;
  //cout << " D(t) = " << _D_t(k_temp,0.01*k_temp) << endl << endl;//*/

  /*double kk = .05; 
  cout << "Dt(" << kk << ") = " <<   _D_t(kk,0.01*kk) << endl;
  kk = .1;
  cout << "Dt(" << kk << ") = " <<   _D_t(kk,0.01*kk) << endl;
  kk = .15;
  cout << "Dt(" << kk << ") = " <<   _D_t(kk,0.01*kk) << endl;
  kk = .2;
  cout << "Dt(" << kk << ") = " <<   _D_t(kk,0.01*kk) << endl; //*/
  //double kk = .85; 
  //cout << "Dt(" << kk << ") = " <<   _D_t(kk,0.01*kk) << endl;

  //cout << " D(t,ps=.5) = " <<  _D_t(1.,.5) << endl;
  //cout << " D(t,ps=.1) = " <<  _D_t(1.,.1) << endl;
  //cout << " D(t,ps=.05) = " <<  _D_t(1.,.05) << endl;
  //cout << " D(t,ps=.01) = " <<  _D_t(1.,.01) << endl;

  //cout << " D(t,ps=.5) = " <<  _D_t(.5,.5) << endl;
  //cout << " D(t,ps=.1) = " <<  _D_t(.5,.1) << endl;
  //cout << " D(t,ps=.05) = " <<  _D_t(.5,.05) << endl;
  //cout << " D(t,ps=.01) = " <<  _D_t(.5,.01) << endl;

  //cout << " D(t,k=0.5,ps=.01) = " << _D_t(.5,.01) << endl;
  //cout << " D(t,k=1.0,ps=.01) = " << _D_t(1.,.01) << endl;
  //cout << " D(t,k=2.0,ps=.01) = " << _D_t(2.,.01) << endl;

  //cout << " 0.5 | D(t) = " << _D_t(0.5,.02) << endl;
  //cout << " 1.0 | D(t) = " << _D_t(1.0,.02) << endl;
  //cout << " 2.0 | D(t) = " << _D_t(2.0,.02) << endl;
  //cout << " 3.0 | D(t) = " << _D_t(3.0,.02) << endl;
  //cout << " 4.0 | D(t) = " << _D_t(4.0,.01) << endl;
  //cout << " 5.0 | D(t) = " << _D_t(5.0,.01) << endl;
  //cout << " 6.0 | D(t) = " << _D_t(6.0,.01) << endl;
  //cout << " 7.0 | D(t) = " << _D_t(7.0,.01) << endl;
  //cout << " 8.0 | D(t) = " << _D_t(8.0,.01) << endl;
  //cout << " 9.0 | D(t) = " << _D_t(9.0,.01) << endl;
  //cout << "10.0 | D(t) = " << _D_t(10.0,.01) << endl;
  //cout << "15.0 | D(t) = " << _D_t(15.0,.01) << endl;

  return 0;
}

/*--------------------------------------------------------------------*/
// assemble masters

struct Rho_V
{
  double lo, nlo;
  int S[3];
  //double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // group factors
  double Nc=1., cF=1.; // QED
  bool calc_nlo = false;

  // notation à la 1310.0164
  Master
    *rho_b, *rho_bb, *rho_d, *rho_db, *rho_g, *rho_gp, *rho_h, *rho_hp, *rho_j;
  double
    _b, _bb, _d, _db, _g, _gp, _h, _hp, _j;

  Rho_V() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    rho_b =  _11020(0,0,S); // 'vector' channel
    rho_bb=  _10120(0,0,S);
    rho_d =  _11010(0,0,S);
    rho_db=  _10110(0,0,S);
    rho_g =  _11011(0,0,S);
    rho_gp=  _11011(1,1,S); // (1,1)
    rho_h =  _11110(0,0,S);
    rho_hp=  _Star( 0,0,S);
    rho_j =  _11111(0,0,S);

  };
  void operator ()() {

    lo = +Nc*K2*psi0(-1,-1,MOT1)*OOFP; // note sign convention
    //lo = Nc*K2*OOFP; // large-K2
    nlo= lo*3.*cF*SQR(OOFP);

    if (calc_nlo) {
    // Quadrature step! --
    double res, err;

    gsl_set_error_handler_off();
    size_t limit = 5e2;

    quad wsp1(limit);
    quad wsp2(limit);

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return (  2.*(rho_h ->integrand)(x,y)
                     +2.*(rho_hp->integrand)(x,y)*(kp/km)
                     -   (rho_j ->integrand)(x,y)*(kp/km)
                   )/SQR(kp);
          } );
      gsl_integration_qag( inner, .0+1e-10,1., tol, 0,
                          limit, 6, wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, .0+1e-10,1., tol*2.,  0,
                        limit, 6, wsp2, &res, &err  );//*/

    // Simpler masters --
    _b = (*rho_b )(k0,k)*K2;
    _bb= (*rho_bb)(k0,k)*K2;
    _d = (*rho_d )(k0,k);
    _db= (*rho_db)(k0,k);
    _g = (*rho_g )(k0,k)*K2*(SQR(k0/k)+3.);
    _gp= (*rho_gp)(k0,k)*K2/SQR(k);
    //_h = (*rho_h )(k0,k)*K2; // these masters are commented b/c I compute them above
    //_hp= (*rho_hp)(k0,k);
    //_j = (*rho_j )(k0,k)*SQR(K2);

    nlo -=
    //8.*Nc*cF*( 2.*(_b-_bb+_d-_db) - 1.5*_g + 2.*(_h+_hp) - _j );
    8.*Nc*cF*( 2.*(_b-_bb+_d-_db) - .5*_g + 2.*_gp + res*CUBE(OOFP)*SQR(kp) );
    } else { nlo = 0.; }

  };
};

struct Rho_00
{
  double lo, nlo;
  int S[3];
  //double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // QCD
  double Nc=1., cF=1.; // QED
  bool calc_nlo = false;

  Master
    *rho_b_0, *rho_bb_0,
    *rho_b_1, *rho_bb_1,
    *rho_b_2, *rho_bb_2,
    *rho_g, *rho_gp,
    *rho_h_0, *rho_h_1,
    *rho_hp,
    *rho_j_0, *rho_j_2;
  double
    _b_1, _bb_1,
    _b_2, _bb_2,
    _g,   _gp,
    _h_0, _h_1,
    _hp,
    _j_0, _j_2;

  Rho_00() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    rho_b_0 =  _11020(0,0,S); // notation: rho_<tag>_<power of ..>
    rho_bb_0=  _10120(0,0,S);
    rho_b_1 =  _11020(1,0,S);
    rho_bb_1=  _10120(1,0,S);
    rho_b_2 =  _11020(2,0,S);
    rho_bb_2=  _10120(2,0,S);
    rho_g   =  _11011(0,0,S);
    rho_gp  =  _11011(1,1,S);
    rho_h_0 =  _11110(0,0,S);
    rho_h_1 =  _11110(0,1,S);
    rho_hp  =  _Star( 0,0,S);
    rho_j_0 =  _11111(0,0,S);
    rho_j_2 =  _11111(2,0,S);

  };
  void operator ()() {

    lo = -2.*Nc*( k0*k0*(psi1(-1,-1,MOT1)-psi2(-1,-1,MOT1))-.25*K2*psi0(-1,-1,MOT1) )*OOFP;
    nlo = -cF*Nc*( k*k*psi0(-1,-1,MOT1) )*CUBE(OOFP);
    //nlo = -cF*Nc*( k*k )*CUBE(OOFP); // large-K2

    if (calc_nlo) {
    // Quadrature step! --
    double res, err;

    gsl_set_error_handler_off();
    size_t limit = 5e2;

    quad wsp1(limit);
    quad wsp2(limit);

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return (  2.*(rho_h_0->integrand)(x,y)
                     +2.*(rho_hp ->integrand)(x,y)*(kp/km)+
                     -8.*(rho_h_1->integrand)(x,y)*(k0*k0)/K2
                     +   (rho_j_0->integrand)(x,y)*(k0*k0+k*k)*(kp/km)/K2
                     -4.*(rho_j_2->integrand)(x,y)*(kp/km)
                   )/SQR(kp);
          } );
      gsl_integration_qag( inner, .0+1e-10,1., tol, 0,
                          limit, 6, wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, .0+1e-10,1., tol*2., 0,
                       limit, 6, wsp2, &res, &err  );//*/

    // Simpler master(s) --
    //_b_0 = (*rho_b_0 )(k0,k)*K2;
    //_bb_0= (*rho_bb_0)(k0,k)*K2;
    //_b_1 = (*rho_b_1 )(k0,k)*k0;
    //_bb_1= (*rho_bb_1)(k0,k)*k0;
    //_b_2 = (*rho_b_2 )(k0,k);
    //_bb_2= (*rho_bb_2)(k0,k);
    _g   = (*rho_g )(k0,k)*( k*k + k0*k0 );
    _gp  = (*rho_gp)(k0,k)*( -4. );
    //_h_0 = (*rho_h_0 )(k0,k)*K2;
    //_h_1 = (*rho_h_1 )(k0,k)*k0;
    //_hp  = (*rho_hp  )(k0,k);
    //_j_0 = (*rho_j_0 )(k0,k)*(k0*k0+k*k)*K2;
    //_j_2 = (*rho_j_2 )(k0,k)*K2;

    nlo -=
    4.*Nc*cF*(// 2.*(_b_0-_bb_0-4.*(_b_1-_bb_1)+4.*(_b_2-_bb_2)) // =0  (GJ: when did I write this???)
                                                                 // (GJ: eq(7.1) of your paper, you loskop!)
              //+ _g + 2.*(_h_0+_hp) - 8.*_h_1 + _j_0 - 4.*_j_2 );
             + _g + _gp + res*CUBE(OOFP)*SQR(kp) );
    } else { nlo=0.; }

  };
};

/*--------------------------------------------------------------------*/
// evaluation

double chi_V(double p0, double p) {
  k0 = fabs(p0); k = fabs(p);
  double P2 = p0*p0 - p*p;
  return 4.*( P2*chi(0,-1,-1) - 1./12. );
}

double chi_00(double p0, double p) {
  k0 = fabs(p0); k = fabs(p);
  double P2 = p0*p0 - p*p;
  return -8.*( SQR(p0)*( chi(2,-1,-1) - chi(1,-1,-1) ) 
                   + .25*( P2*chi(0,-1,-1) - 1./12. ) );
}

double chi_T(double p0, double p) {
  double P2 = p0*p0 - p*p;
  return -.5*( chi_V(p0,p) - chi_00(p0,p)*P2/(SQR(p)) );
}

double chi_L(double p0, double p) {
  double P2 = p0*p0 - p*p;
  return - chi_00(p0,p)*P2/(SQR(p));
}

// HTL limits:

double rho_T_HTL(double p0, double p) {
  double P2 = p0*p0 - p*p;
  double p3 = p*p*p;
  if (fabs(p0)<p) { return + M_PI*p0*P2/(12.*p3); }
  else { return 0.; }
}

double rho_L_HTL(double p0, double p) {
  double P2 = p0*p0 - p*p;
  double p3 = p*p*p;
  if (fabs(p0)<p) { return - M_PI*p0*P2/(6.*p3); }
  else { return 0.; }
}

double chi_T_HTL(double p0, double p) {
  double P2 = p0*p0 - p*p;
  double LG = log(fabs((p+p0)/(p-p0)));
  return ( p0*p0 - .5*(p0*P2/p)*LG )/(6.*p*p)
    ; 
}

double chi_L_HTL(double p0, double p) {
  double P2 = p0*p0 - p*p;
  double LG = log(fabs((p+p0)/(p-p0)));
  return - ( 1. - .5*(p0/p)*LG )*P2/(3.*p*p)
    ;
}

/*--------------------------------------------------------------------*/

#include "gauss.h"

double alpha_em = 1./137.0359991;

/*--------------------------------------------------------------------*/

//#include "gauss.h"
//
double _A_integrand(double kk, double p_m, double p_p) {
  // LO part
  double p0 = p_p + p_m;
  double p  = p_p - p_m;
  double P2 = 4.*p_p*p_m;

  Rho_V rV;
  Rho_00 r00;
  CHEM_POT=false;
  MOT1 = 0.; MOT2 = 0.;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;

  double sign_k0 = 1.;
  if (p0 < 0.) { sign_k0 = -1.; }

  k0 = fabs(p0); k = fabs(p); // NB: k != kk

  double tkmp = SQR( (2*kk - p0)/p );

  rV(); r00();

  double V_lo  = rV.lo;
  double L_lo  = - r00.lo*P2/p/p;
  double T_lo  = -.5*(- V_lo + L_lo);

  double res_lo  = .5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );

  res_lo *=  (rV.lo)*(1.+tkmp) -  r00.lo*P2*(1.-3.*tkmp)/p/p;

  return -sign_k0*res_lo;

}

double _B_integrand(double kk, double p_m, double p_p) {
  // LO part
  double p0 = p_p + p_m;
  double p  = p_p - p_m;
  double P2 = 4.*p_p*p_m;

  Rho_V rV;
  Rho_00 r00;
  CHEM_POT=false;
  MOT1 = 0.; MOT2 = 0.;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;

  double sign_k0 = 1.;
  if (p0 < 0.) { sign_k0 = -1.; }

  k0 = fabs(p0); k = fabs(p); // NB: k != kk

  double tkmp = SQR( (2*kk - p0)/p );

  rV(); r00();

  double V_lo  = rV.lo;
  double L_lo  = - r00.lo*P2/p/p;
  double T_lo  = -.5*(- V_lo + L_lo);

  double res_lo  = .5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );

  res_lo *=  (rV.lo)*(1.+tkmp) -  r00.lo*P2*(1.-3.*tkmp)/p/p;

  return -sign_k0*res_lo*log(fabs(P2));

}


double _C_integrand(double kk, double p_m, double p_p) {
  // LO part
  double p0 = p_p + p_m;
  double p  = p_p - p_m;
  double P2 = 4.*p_p*p_m;

  Rho_V rV;
  Rho_00 r00;
  CHEM_POT=false;
  MOT1 = 0.; MOT2 = 0.;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;

  rV.calc_nlo = true;
  r00.calc_nlo = true;

  double sign_k0 = 1.;
  if (p0 < 0.) { sign_k0 = -1.; }

  k0 = fabs(p0); k = fabs(p); // NB: k != kk

  double tkmp = SQR( (2*kk - p0)/p );

  rV(); r00();

  double V_nlo = rV.nlo;
  double L_nlo = - r00.nlo*P2/p/p;
  double T_nlo = -.5*(- V_nlo + L_nlo);

  double res_nlo = .5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );

  res_nlo *= (rV.nlo)*(1.+tkmp) - r00.nlo*P2*(1.-3.*tkmp)/p/p;

  return -sign_k0*res_nlo;

}

double _D_integrand(double kk, double p_m, double p_p) {
  // LO part
  double p0 = p_p + p_m;
  double p  = p_p - p_m;
  double P2 = 4.*p_p*p_m;

  Rho_V rV;
  Rho_00 r00;
  CHEM_POT=false;
  MOT1 = 0.; MOT2 = 0.;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;

  double sign_k0 = 1.;
  if (p0 < 0.) { sign_k0 = -1.; }

  k0 = fabs(p0); k = fabs(p); // NB: k != kk

  double tkmp = SQR( (2*kk - p0)/p );

  rV(); r00();

  double V_lo  = rV.lo;
  double V_nlo = rV.nlo;
  double L_lo  = - r00.lo*P2/p/p;
  double L_nlo = - r00.nlo*P2/p/p;
  double T_lo  = -.5*(- V_lo + L_lo);
  double T_nlo = -.5*(- V_nlo + L_nlo);
  double T_chi = chi_T(p0,p);
  double L_chi = chi_L(p0,p);

  double counter = p*( 1.+f(kk-p0,-1)+f(p0,+1) );

  counter  *= 2.*( T_lo*T_chi + L_lo*L_chi + (T_lo*T_chi-L_lo*L_chi)*tkmp );

  return -sign_k0*counter;

}



void QED_integrand_(double kk, double p_m, double p_p, double *A, double *B, double *C, double *D) {

  double p0 = p_p + p_m;
  double p  = p_p - p_m;
  double P2 = 4.*p_p*p_m;

  Rho_V rV;
  Rho_00 r00;
  CHEM_POT=false;
  MOT1 = 0.; MOT2 = 0.;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;

  double sign_k0 = 1.;
  if (p0 < 0.) { sign_k0 = -1.; }

  k0 = fabs(p0); k = fabs(p); // NB: k != kk

  double tkmp = SQR( (2*kk - p0)/p );

  rV(); r00();

  double V_lo  = rV.lo;
  double V_nlo = rV.nlo;
  double L_lo  = - r00.lo*P2/p/p;
  double L_nlo = - r00.nlo*P2/p/p;
  double T_lo  = -.5*(- V_lo + L_lo);
  double T_nlo = -.5*(- V_nlo + L_nlo);
  double T_chi = chi_T(p0,p);
  double L_chi = chi_L(p0,p);

  double res_lo  = .5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );
  double res_nlo = .5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );
  double counter = p*( 1.+f(kk-p0,-1)+f(p0,+1) );

  res_lo *=  (rV.lo)*(1.+tkmp) -  r00.lo*P2*(1.-3.*tkmp)/p/p;
  res_nlo *= (rV.nlo)*(1.+tkmp) - r00.nlo*P2*(1.-3.*tkmp)/p/p;
  counter  *= 2.*( T_lo*T_chi + L_lo*L_chi + (T_lo*T_chi-L_lo*L_chi)*tkmp );

  *A = sign_k0*res_lo;
  *B = sign_k0*res_lo*log(fabs(P2));
  *C = sign_k0*res_nlo;
  *D = sign_k0*counter;

}

double integrator(double x1, double x2, double y1, double y2, double (*func)(double,double,double),void (*params)) {
  double x, y, res = 0., temp;
  double k      = ((double *)params)[0];
  double p_star = ((double *)params)[1];
  double ai, ti, aj, tj;
  for (int i=0;i<64;i++) { 
    ai  = .5*G64pt[i][0]; ti = .5*(G64pt[i][1]+1.);
    ai *= x2-x1;
    x   = x1*(1.-ti) + x2*ti;
    for (int j=0;j<64;j++) { 
      aj  = .5*G64pt[j][0]; tj = .5*(G64pt[j][1]+1.);
      aj *= y2-fmax(y1,p_star+x);
      y   = fmax(y1,p_star+x)*(1.-tj) + y2*tj;
      temp = func(k,x,y);
      res += ai*aj*temp;
    }
  }
  return res;
}

double region_s(double k,double (*func)(double,double,double)) {
  double k_in[2] = { k, 0. };
  double res;
  res = -integrator(.0,k,k,27.,func,k_in);
  return res/CUBE(k);
}

double region_t(double k,double (*func)(double,double,double)) {
  double k_in[2] = { k, 0. };
  double res;
  double e2 = (4.*M_PI*alpha_em);
  res = integrator(-14.013-k*.5,0.,0.,k,func,k_in);
  return res/CUBE(k) ;
}

double _D_t(double k, double p_star) {
  //double p_star = .5;
  double k_in[2] = { k, p_star };
  //double res;
  double e = sqrt(4.*M_PI*alpha_em);
  //double p_star = .5;
    // Quadrature step! --
    double res, err;

    gsl_set_error_handler_off();
    size_t limit = 5e1;

    quad wsp1(limit);
    quad wsp2(limit);

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return _D_integrand(k,x,y);
          } );
      gsl_integration_qag( inner, fmax(1e-9,p_star+x), k, tol, 0,
                          limit, 6, wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, -11.-.4*k, 1e-9, 0, tol*2.,
                        limit, 6, wsp2, &res, &err  );//*/

  //res  = integrator(-k*.5,0.,0.,k,_D_integrand,k_in);
  //res += integrator(-1.,-k*.5,0.,k,_D_integrand,k_in);
  //res += -k*k*(2.*M_PI/9.)*( log(p_star/e) + .2181 );
  res += -k*k*(2.*M_PI/9.)*( log(p_star/e) + .554521 );
  //res += -k*k*(2.*M_PI/9.)*( log(p_star) + .56 );
  return res/CUBE(k);
}

double _A_t(double k) {
  double e = sqrt(4.*M_PI*alpha_em);
    double res, err;

    gsl_set_error_handler_off();
    size_t limit = 5e1;

    quad wsp1(limit);
    quad wsp2(limit);

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return _A_integrand(k,x,y);
          } );
      gsl_integration_qag( inner, 1e-9, k, tol, 0,
                          limit, 6, wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, -17.-.4*k, 1e-9, 0, tol*2.,
                        limit, 6, wsp2, &res, &err  );//*/

  return res/CUBE(k);
}

double _B_t(double k) {
  double e = sqrt(4.*M_PI*alpha_em);
    double res, err;

    gsl_set_error_handler_off();
    size_t limit = 5e1;

    quad wsp1(limit);
    quad wsp2(limit);

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return _B_integrand(k,x,y);
          } );
      gsl_integration_qag( inner, 1e-9, k, tol, 0,
                          limit, 6, wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, -17.-.4*k, 1e-9, 0, tol*2.,
                        limit, 6, wsp2, &res, &err  );//*/

  return res/CUBE(k);
}





int Print_nu_rate(double k_min,double k_max,int N_k,string fname,double mu=0.) {
  double res, s;
  s=(k_max-k_min)/((double)N_k-1.);

  // filename
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  fout.open(fname);
  fout << "# Columns: k/T, A(s), A(t), B(s), B(t), C(s), C(t), D(s), D(t) " 
       << endl
       << "# ( where A, B, C, D are defined in 2312.07015 )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  double kk=k_min;

  double res_s_A=0., res_s_B=0., res_s_C=0., res_s_D=0.,
         res_t_A=0., res_t_B=0., res_t_C=0., res_t_D=0.;

  for (int i=0;i<N_k;i++) { 

    res_s_A = region_s(kk,_A_integrand);
    res_t_A = region_t(kk,_A_integrand);

    res_s_B = region_s(kk,_B_integrand);
    res_t_B = region_t(kk,_B_integrand);

    res_s_C = region_s(kk,_C_integrand);
    res_t_C = region_t(kk,_C_integrand);

    res_s_D = region_s(kk,_D_integrand);
    res_t_D = _D_t(kk,0.01*kk); // p* = k/100 [soft/hard cut off]

    percentage = (kk-k_min)/(k_max-k_min);
    fout << scientific << kk
         <<     "    " << res_s_A
         <<     "    " << res_t_A
         <<     "    " << res_s_B
         <<     "    " << res_t_B
         <<     "    " << res_s_C
         <<     "    " << res_t_C
         <<     "    " << res_s_D
         <<     "    " << res_t_D
         << endl;

    kk+=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

