#include "core.hh"
#include "quad.hh"
#include "map.hh"
#include "timer.hh"
#include <fstream>
#include <string>

double k0,k;
double mu_q = 0.;
bool  CHEM_POT = false;
double MOT1, MOT2, MOT3, MOT4, MOT5;// mu over T
double tol = 1e-1; // absolute tolerance

using namespace std;

ofstream fout;
ifstream fin;
int elapsed; float percentage;
void D(double,double);

// QED
int Print_QED(double,double,int,string,double);
void QED_integrand_(double kk, double p_m, double p_p, double *_lo, double *_nlo);

// lattice
int ReadIn(string,string,double,double);

// heavy ion rates:
int hydro_table_integrated(string);
int hydro_table_unintegrated(string);
int hydro_table_integrated_T_L(string,string);

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

  Print_QED(5.,6.,10,"QED_rho_7",0.);
  // Return: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
  //D(k0,k);

}

/*--------------------------------------------------------------------*/
// assemble masters

struct Rho_V
{
  double lo, nlo;
  int S[3];
  //double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // group factors
  double Nc=1., cF=1.; // QED

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

  };
};

struct Rho_00
{
  double lo, nlo;
  int S[3];
  //double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // QCD
  double Nc=1., cF=1.; // QED

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

  };
};

/*--------------------------------------------------------------------*/
// evaluation

void D(double k0_curr, double k_curr) {
  k0=k0_curr; k=k_curr;

  Rho_V   rV; // assign ptrs
  Rho_00 r00;

  rV(); r00();


  cout << scientific << k0        // k0/T
       <<     "    " << rV.lo     // leading-order: rho_V ,
       <<     "    " << r00.lo    //                rho_00
       <<     "    " << rV.nlo    // next-to-LO   : rho_V ,
       <<     "    " << r00.nlo   //                rho_00
       << endl;
}

/*--------------------------------------------------------------------*/

#include "gauss.h"

// some physical units, at last:

double    hbarc = .19732698041522198110 ; // GeV.fm
double alpha_em = 1./137.0359991;
double      m_e = .0005109895, // GeV
            m_m = .10565837;   // GeV

/*--------------------------------------------------------------------*/

//#include "gauss.h"

void QED_integrand_(double kk, double p_m, double p_p, double *_lo, double *_nlo) {

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

  cout << " p0 = " << p0 << " , p = " << p << endl;
  cout << " rV = " << sign_k0*rV.lo << " , r00 = " << -sign_k0*r00.lo << endl;

  double V_lo  = rV.lo;
  double V_nlo = rV.nlo;
  double T_lo  = - r00.lo*P2/p/p;
  double T_nlo = - r00.nlo*P2/p/p;
  double L_lo  = V_lo - 2.*T_lo;
  double L_nlo = V_nlo - 2.*T_nlo;

  double res_lo  = -.5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );
  double res_nlo = -.5*p*P2*( 1.+f(kk-p0,-1)+f(p0,+1) );

//  res_lo *=  (rV.lo)*(1.+tkmp) -  r00.lo*P2*(1.-3.*tkmp)/p/p;
//  res_nlo *= (rV.nlo)*(1.+tkmp) - r00.nlo*P2*(1.-3.*tkmp)/p/p;
  res_lo  *= 2.*( T_lo + L_lo + (T_lo-L_lo)*tkmp );
  res_nlo *= 2.*( T_nlo+ L_nlo+ (T_nlo-L_nlo)*tkmp );

  *_lo  = sign_k0*res_lo;
  *_nlo = sign_k0*res_nlo;

}

int Print_QED(double k_min,double k_max,int N_k,string fname,double mu=0.) { // k = y.T, k0 = x.y.T
  //int N_k;
  int N_pm, N_pp;
  double res, s;// k_min, k_max;
  //k=y;
  s=(k_max-k_min)/((double)N_k-1.);

  // filename
  //string fname = "QED_rho11";
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  fout.open(fname);
  fout << 
  "# Columns: k/T, Gamma_lo/GF2, Gamma_nlo/(e2*GF2)" 
       << endl
       << "# ( k=" << k << " )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  //int N_pm=10; 
  //int N_pp=10; 
  double pm, pp;

  // first: t-channel
  double pm_min;
  double pm_max;
  double pp_min;
  double pp_max; // depends on kk
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  double s_m = (pm_max-pm_min)/((double)N_pm-1.);
  double s_p = (pp_max-pp_min)/((double)N_pp-1.);
  //s = 1e-1;
  //k0+=s;
  double kk=k_min;

  double res_lo=0., res_nlo=0., _lo,_nlo;

  double ai_m, ti_m, ai_p, ti_p;
  for (int i=0;i<N_k;i++) { 

    res_lo = 0.;
    res_nlo = 0.;

    // do t-channel integral
    pm_min=-10.;
    pm_max=-.01;
    pp_min=.01;
    pp_max=kk; // fix k
    for (int i_m=0;i_m<32;i_m++) { 
      ai_m  = .5*G32pt[i_m][0]; ti_m = .5*(G32pt[i_m][1]+1.);
      ai_m *= pm_max-pm_min;
      pm    = pm_min*(1.-ti_m) + pm_max*ti_m;
      for (int i_p=0;i_p<32;i_p++) { 
        ai_p  = .5*G32pt[i_p][0]; ti_p = .5*(G32pt[i_p][1]+1.);
        ai_p *= pp_max-pp_min;
        pp    = pp_min*(1.-ti_p) + pp_max*ti_p;
        cout << " pp = " << pp << " , pm = " << pm << endl;
        QED_integrand_(kk,pm,pp,&_lo,&_nlo);
        cout << " ... point complete ..." << endl;
        res_lo  += ai_p*ai_m*_lo;
        res_nlo += ai_p*ai_m*_nlo;
      }
    }

    // do s-channel integral
    pm_min=.01;
    pm_max=kk; // fix k
    pp_min=kk;
    pp_max=14.;
    for (int i_m=0;i_m<32;i_m++) { 
      ai_m  = .5*G32pt[i_m][0]; ti_m = .5*(G32pt[i_m][1]+1.);
      ai_m *= pm_max-pm_min;
      pm    = pm_min*(1.-ti_m) + pm_max*ti_m;
      for (int i_p=0;i_p<32;i_p++) { 
        ai_p  = .5*G32pt[i_p][0]; ti_p = .5*(G32pt[i_p][1]+1.);
        ai_p *= pp_max-pp_min;
        pp    = pp_min*(1.-ti_p) + pp_max*ti_p;
        cout << " pp = " << pp << " , pm = " << pm << endl;
        QED_integrand_(kk,pm,pp,&_lo,&_nlo);
        cout << " ... point complete ..." << endl;
        res_lo  -= ai_p*ai_m*_lo;
        res_nlo -= ai_p*ai_m*_nlo;
      }
    }

    percentage = (kk-k_min)/(k_max-k_min);
    fout << scientific << kk
         <<     "    " << res_lo
         <<     "    " << res_nlo
         << endl;

    kk+=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}
