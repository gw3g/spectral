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
int Print_D(double); int elapsed; float percentage;
int Print_D(double,double);
int Print_HTL(double,double);
int Print2_HTL(double,double);
int Print_k2ave(double); 
void D(double,double);

// lattice
int ReadIn(string,string,double,double);

// heavy ion rates:
int hydro_table_integrated(string);
int hydro_table_unintegrated(string);
int hydro_table_integrated_T_L(string,string);

int main(int argc, char *argv[]) {

  // Enter: k0, k
  if (argc>1) { k0=atof(argv[1]); k=atof(argv[2]); }

  // Return: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
  //D(k0,k);

  /* for lattice comparisons: */
  //nf=0
  //Print_D(3.*2.*M_PI/3.,0.);     // T=1.1Tc
  //Print_D(1.*2.*M_PI/3.,3./3.);     // T=1.1Tc
  //Print_D(1.*2.*M_PI/3.,3./3.);     // T=1.1Tc
  //Print_k2ave(5.);     // T=1.1Tc
  //Print_D(3.*7.*M_PI/12.);    // T=1.3Tc
  // nf=2
  //Print_D(sqrt(1.)*M_PI/2.); // T=1.2Tc
  //Print_D(1.5*M_PI);         // T=1.2Tc
  //Print_D(.5*M_PI*sqrt(14)); // T=1.2Tc
  //
  //
  char  name_in[100];
  char  name_out[100];
  char *tag;
  double MUq;

  /* fixed alpha */

  tag = ".R1(5)";
  int Nf=3; k = 2.*M_PI; //double t=1.8; 
  MUq = 12./3.;
  sprintf(name_in,"fixed_coupling.dat");
  sprintf(name_out,"out_for_plot/NLO_rho_{k=%.2f,mu=%.2f}.dat",k,MUq);
  ReadIn(string(name_in),string(name_out),k,MUq);//*/

  /* lattice table */
  /*

  tag = ".R1(5)";
  int Nf=3; k = 1*2.*M_PI/3.; double t=1.8; 
  MUq = 3./3.;
  sprintf(name_in ,"coupling_nf%d_{k=%.2f,t=%.2f,mu=%.2f}%s.dat",Nf,k,t,MUq,tag);
  sprintf(name_out,"out/NLO_rho_{k=%.2f,mu=%.2f}.dat",Nf,k,MUq);
  ReadIn(string(name_in),string(name_out),k,MUq);*/


  /* for hydro run: */
  //hydro_table_T_L();
  //hydro_table_integrated_T_L("meshK_alpha2.dat",
  //                           "meshK_NLO_kT_integrated_T_L_2.dat");

  /* Others: */
  //Print_D(.02,2.);
  //Print_D(.25,10.);
  //Print_D(.25,1.);
  //Print_D(.25,2.);

  //Print_HTL(1e-1,0.);
  //Print_HTL(1e-2,2.);
  //Print_HTL(1e-3,0.);

  //Print2_HTL(.3,0.01);
  //Print2_HTL(.4,0.01);

  //Print_D(.3);
  //Print_D(1.5);
  //Print_D(3.);
  //Print_D(6.);
  //Print_D(9.);
}

/*--------------------------------------------------------------------*/
// assemble masters

struct Rho_V
{
  double lo, nlo;
  int S[3];
  double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // group factors
  //double Nc=1., cF=1.; // QED

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
  double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // QCD
  //double Nc=1., cF=1.; // QED

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

int Print_D(double k_curr,double mu=0.) {
  int N_k0;
  double res, s, k0_min, k0_max;
  k=k_curr;

  CHEM_POT=true;
  MOT1 = mu; MOT2 = -mu;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
  // filename
  char k_name[20];
  char mu_name[20];
  sprintf(k_name,"{k=%.2f}",k);
  sprintf(mu_name,"{mu=%.2f}",mu);
  string fname = "NLO_rho_"
               + string(k_name);
  if (CHEM_POT) {fname = fname+"."+string(mu_name); }
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV, rrV;
  Rho_00 r00, rr00;
  fout.open(fname);
  fout << 
  "# Columns: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
       << endl
       << "# ( k=" << k << " )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_k0=50; 

  k0_min=1e-4;
  k0_max=.5;
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  s=(k0_max-k0_min)/((double)N_k0-1.);
  //s = 1e-1;
  //k0+=s;
  k0=k0_min;

  int i=0;
  //for (int i=0;i<N_k0;i++) { 
  while (k0<k0_max) {
    //percentage=(float)i/((float)N_k0);
    percentage = (k0-k0_min)/(k0_max-k0_min);

    MOT1 = mu; MOT2 = -mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV(); r00();
    MOT1 = -mu; MOT2 = mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rrV(); rr00();

    fout << scientific << k0                    // k0/T
         <<     "    " << .5*(rV.lo+rrV.lo)     // leading-order: rho_V ,
         <<     "    " << .5*(r00.lo+rr00.lo)   //                rho_00
         <<     "    " << .5*(rV.nlo+rrV.nlo)   // next-to-LO   : rho_V ,
         <<     "    " << .5*(r00.nlo+rr00.nlo) //                rho_00
         << endl;

    k0+=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

int Print_HTL(double y,double mu=0.) { // k = y.T, k0 = x.y.T
  int N_x;
  double res, s, x_min, x_max, x;
  k=y;

  CHEM_POT=true;
  MOT1 = mu; MOT2 = -mu;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
  // filename
  char k_name[20];
  char mu_name[20];
  sprintf(k_name,"{y=%.3f}",y);
  sprintf(mu_name,"{mu=%.2f}",mu);
  string fname = "NLO_rho_"
               + string(k_name);
  if (CHEM_POT) {fname = fname+"."+string(mu_name); }
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV, rrV;
  Rho_00 r00, rr00;
  fout.open(fname);
  fout << 
  "# Columns: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
       << endl
       << "# ( k=" << k << " )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_x=50; 

  x_min=.1;
  x_max=2.;
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  s=(x_max-x_min)/((double)N_x-1.);
  //s = 1e-1;
  //k0+=s;
  x=x_min;

  int i=0;
  //for (int i=0;i<N_k0;i++) { 
  while (x<x_max) {
    //percentage=(float)i/((float)N_k0);
    percentage = (x-x_min)/(x_max-x_min);
    k0 = x*y;

    MOT1 = mu; MOT2 = -mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV(); r00();
    MOT1 = -mu; MOT2 = mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rrV(); rr00();

    fout << scientific << x                     // k0/k
         <<     "    " << .5*(rV.lo+rrV.lo)     // leading-order: rho_V ,
         <<     "    " << .5*(r00.lo+rr00.lo)   //                rho_00
         <<     "    " << .5*(rV.nlo+rrV.nlo)   // next-to-LO   : rho_V ,
         <<     "    " << .5*(r00.nlo+rr00.nlo) //                rho_00
         << endl;

    x+=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}



int Print2_HTL(double x,double y) { // k = y.T, k0 = x.y.T
  int N_mu;
  double res, s, mu_min, mu_max, mu;
  k=y;
  k0=x*y;

  CHEM_POT=true;
  //MOT1 = mu; MOT2 = -mu;
  //MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
  // filename
  char k_name[20];
  //char mu_name[20];
  sprintf(k_name,"{x=%.3f,y=%.3f}",x,y);
  //sprintf(mu_name,"{mu=%.2f}",mu);
  string fname = "NLO_rho_"
               + string(k_name);
  if (CHEM_POT) {fname = fname; }
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV, rrV;
  Rho_00 r00, rr00;
  fout.open(fname);
  fout << 
  "# Columns: mu/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
       << endl
       << "# ( k=" << mu << " )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_mu=60; 

  mu_min=.0;
  mu_max=30.;
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  s=(mu_max-mu_min)/((double)N_mu-1.);
  //s = 1e-1;
  //k0+=s;
  mu=mu_min;

  int i=0;
  //for (int i=0;i<N_k0;i++) { 
  while (mu<mu_max) {
    //percentage=(float)i/((float)N_k0);
    percentage = (mu-mu_min)/(mu_max-mu_min);

    MOT1 = mu; MOT2 = -mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV(); r00();
    MOT1 = -mu; MOT2 = mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rrV(); rr00();

    fout << scientific << mu                    // mu/T
         <<     "    " << .5*(rV.lo+rrV.lo)     // leading-order: rho_V ,
         <<     "    " << .5*(r00.lo+rr00.lo)   //                rho_00
         <<     "    " << .5*(rV.nlo+rrV.nlo)   // next-to-LO   : rho_V ,
         <<     "    " << .5*(r00.nlo+rr00.nlo) //                rho_00
         << endl;

    mu+=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

/*--------------------------------------------------------------------*/

#include <gsl/gsl_sf_bessel.h>
double k2av(double M) {
  // worth noting: k0=3*k for k0~26.
  return 3.*M*gsl_sf_bessel_Kn(3,M)/gsl_sf_bessel_Kn(2,M);
}

int Print_k2ave(double mu=0.) {
  int N_M;
  double res, s, M, M_min, M_max;

  CHEM_POT=true;
  MOT1 = mu; MOT2 = -mu;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
  // filename
  char k_name[20];
  char mu_name[20];
  sprintf(k_name,"k2ave");
  sprintf(mu_name,"{mu=%.2f}",mu);
  string fname = "NLO_rho_"
               + string(k_name);
  if (CHEM_POT) {fname = fname+"."+string(mu_name); }
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV, rrV;
  Rho_00 r00, rr00;
  fout.open(fname);
  cout << 
  "# Columns: M/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
       << endl
       << "# ( mu=" << mu << " )" << endl;

  //signal( SIGALRM, sigalrm_handler );
  //elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_M=20; 

  M_min=1e-1;
  M_max=1e2;
  // don't change anything after that.

  s=pow(M_max/M_min,1./(N_M-1));
  //s=(k0_max-k0_min)/((double)N_k0-1.);
  //s = 1e-1;
  //M*=s;
  M=M_min;

  for (int i=0;i<N_M;i++) { 
  //while (M<M_max) {
    percentage=(float)i/((float)N_M);
    //percentage = (M-M_min)/(M_max-M_min);
    k=sqrt( fabs(k2av(M)) );
    k0=sqrt(M*M+k*k);

    cout << "percentage = " << 100*percentage << endl;
    MOT1 = mu; MOT2 = -mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV(); r00();
    MOT1 = -mu; MOT2 = mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rrV(); rr00();

    fout << scientific << M                     // M/T
         <<     "    " << .5*(rV.lo+rrV.lo)     // leading-order: rho_V ,
         <<     "    " << .5*(r00.lo+rr00.lo)   //                rho_00
         <<     "    " << .5*(rV.nlo+rrV.nlo)   // next-to-LO   : rho_V ,
         <<     "    " << .5*(r00.nlo+rr00.nlo) //                rho_00
        << endl;

    M*=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

/*--------------------------------------------------------------------*/

int ReadIn(string _IN, string _OUT, double k_curr, double mu=0.) {
  double k0_curr;
  double mubar, aM, aL, aH; // <-- not needed
  k=k_curr;

  CHEM_POT=true;
  MOT1 = mu; MOT2 = -mu;
  MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
  // filename
  fin.open(_IN);
  cout << " Opening input file [" << _IN << "]" << endl;
  fin.ignore(64,'\n');
  fin.ignore(64,'\n');
  fin.ignore(64,'\n');

  fout.open(_OUT);
  /*char opt = 'O';
  if (fout) {
    fout.ignore(64,'\n');
    fout.ignore(64,'\n');
    fout.ignore(64,'\n');
    cout << "\n ... output file exists: overwrite (O) or update (U)? " << endl;
    cin  >> opt;
  }//*/

  /*if ((opt=='U')||(opt=='u')) { // UPDATE
    printf(" Updating.\n");
    while (!fin.eof()) {
      fin  >> k0 >> mubar >> aM >> aL >> aH; fin.ignore(64,'\n');
      fout >> k0_curr >> mubar >> aM >> aL >> aH; fout.ignore(64,'\n');
      if (k0==k0_curr) {
        printf(" done: k0=%f\n",k0_curr);
      }
      else {
        printf("Mismatch! Cannot update, exiting ... \n");
        return;
      }
    }
  }//*/

  //if ((opt=='O')||(opt=='o')) { // OVERWRITE
    k0 = -1.;
    fout << 
    "# Columns: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
         << endl
         << "# ( k=" << k << " )" << endl;

  //}

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV, rrV;
  Rho_00 r00, rr00;

  //signal( SIGALRM, sigalrm_handler );
  //elapsed=0; alarm(1);

  while (k0<20.) {
    fin >> k0 >> mubar >> aM >> aL >> aH; fin.ignore(64,'\n');
    cout << "k0 = " << k0 << endl;

    MOT1 = mu; MOT2 = -mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV(); r00();
    MOT1 = -mu; MOT2 = mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rrV(); rr00();

    fout << scientific << k0                    // k0/T
         <<     "    " << .5*(rV.lo+rrV.lo)     // leading-order: rho_V ,
         <<     "    " << .5*(r00.lo+rr00.lo)   //                rho_00
         <<     "    " << .5*(rV.nlo+rrV.nlo)   // next-to-LO   : rho_V ,
         <<     "    " << .5*(r00.nlo+rr00.nlo) //                rho_00
         << endl;

    //k0+=s; 
  }
  cout << endl << ":: Saved to file [" << _OUT << "]" << endl;
  fin.close();
  fout.close();

  return 0;
}


/*--------------------------------------------------------------------*/

#include "gauss.h"

// some physical units, at last:

double    hbarc = .19732698041522198110 ; // GeV.fm
double alpha_em = 1./137.0359991;
double      m_e = .0005109895, // GeV
            m_m = .10565837;   // GeV

double B_(double x) {
  if (4.*x>1.) return 0.;
  else return (1.+2.*x)*sqrt(1.-4.*x);
}

int hydro_table_integrated(string mesh_name) {
  int N_kT;
  double res_e, res_m, s, kT_min=0., kT_max=15.;
  double y = 0.;
  double alpha, B_e, B_mu;
  double ai, ti;

  string fname = "mesh_GJ_kT_integrated.dat";
  fout.open(fname);
  fout << "# Columns: T/GeV, M/GeV, rate" << endl;

  Rho_V rV;
  Rho_00 r00;

  double T,M;
  fin.open(mesh_name);
  //for (int i=0;i<32000;i++) 
  fin.ignore(64,'\n');

  while (!fin.eof()) {
    fin >> T >> M >> alpha;
    fin.ignore(64,'\n');

    cout << "T=" << T << ", M=" << M << ", alpha=" << alpha << endl;
    B_e  = B_(m_e*m_e/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.); // units:
    B_mu = B_(m_m*m_m/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.); // GeV^-4.fm^-4

    M = M/T;

    res_e = 0.;
    res_m = 0.;
    for (int i=0;i<32;i++) {
      ai  = .5*G32pt[i][0]; ti = .5*(G32pt[i][1]+1.);
      ai *= kT_max-kT_min;
      k   = kT_min*(1.-ti) + kT_max*ti;// = k/T
      k0  = sqrt( M*M+k*k );
      rV();
      res_e += ai*k*( rV.lo + alpha*4.*M_PI*rV.nlo )*exp(-k0)/(1.-exp(-k0));
      res_m += ai*k*( rV.lo + alpha*4.*M_PI*rV.nlo )*exp(-k0)/(1.-exp(-k0));
    }
    res_e *= 2.*M_PI*pow(T,3.)*B_e/(3.*pow(M_PI,3.)*M);
    res_m *= 2.*M_PI*pow(T,3.)*B_mu/(3.*pow(M_PI,3.)*M);

    fout << scientific << T
         <<     "    " << M*T
         <<     "    " << res_e
         <<     "    " << res_m
         << endl;

    //cout << "temp/GeV = " << T << ", M/GeV=" << M << ", k/GeV=" << k << endl;
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();//*/

  return 0;
}


int hydro_table_integrated_T_L(string mesh_name, string fname) {
  int N_kT;
  double res_T_e, res_T_m, 
         res_L_e, res_L_m,
         s, kT_min=0., kT_max=15.;
  double y = 0.;
  double alpha, B_e, B_mu, prefactor;
  double ai, ti;

  //string fname = "mesh_NLO_kT_integrated_T_L.dat";
  fout.open(fname);
  cout << " Opening input file [" << fname << "]\n\n";
  //fout << "# Columns: T/GeV, M/GeV, rate_T(e+e-), rate_T(mu+mu-), rate_L(e+e-), rate_L(mu+mu-)" << endl;
  fout << "# Columns: T/GeV, muB/GeV, M/GeV, rate_T(e+e-), rate_T(mu+mu-), rate_L(e+e-), rate_L(mu+mu-)" << endl;

  Rho_V rV;
  Rho_00 r00;

  CHEM_POT=true;
  Rho_V rrV;
  Rho_00 rr00;

  double T,muB,M, M2, k2, g2, emk0;
  fin.open(mesh_name);
  //for (int i=0;i<32000;i++) 
  fin.ignore(64,'\n');

  while (!fin.eof()) {
    fin >> T >> muB >> M >> alpha;
    fin.ignore(64,'\n');

    cout << "  T = " << T << " ,  mu_B = " << muB << " ,  M = " << M << endl;
    B_e  = B_(m_e*m_e/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.); // units:
    B_mu = B_(m_m*m_m/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.); // GeV^-4.fm^-4

    M = M/T;

    res_T_e = 0.;
    res_T_m = 0.;
    res_L_e = 0.;
    res_L_m = 0.;

    for (int i=0;i<16;i++) {
      ai  = .5*G16pt[i][0]; ti = .5*(G16pt[i][1]+1.);
      ai *= kT_max-kT_min;
      k   = kT_min*(1.-ti) + kT_max*ti;// = k/T
      k2  = k*k; M2 = M*M;
      k0  = sqrt( M2+k2 );
      g2  = alpha*4.*M_PI;
      emk0= exp(-k0);

      MOT1 = muB/(3.*T); MOT2 = -muB/(3.*T);
      MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
      rV();
      r00();

      MOT1 = -muB/(3.*T); MOT2 = muB/(3.*T);
      MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
      rrV();
      rr00();


      res_T_e += ai*k*.25*( rV.lo + r00.lo*(M2/k2) + g2*( rV.nlo + r00.nlo*(M2/k2) ) )*emk0/(1.-emk0);
      res_T_e += ai*k*.25*(rrV.lo +rr00.lo*(M2/k2) + g2*(rrV.nlo +rr00.nlo*(M2/k2) ) )*emk0/(1.-emk0);

      res_T_m += ai*k*.25*( rV.lo + r00.lo*(M2/k2) + g2*( rV.nlo + r00.nlo*(M2/k2) ) )*emk0/(1.-emk0);
      res_T_m += ai*k*.25*(rrV.lo +rr00.lo*(M2/k2) + g2*(rrV.nlo +rr00.nlo*(M2/k2) ) )*emk0/(1.-emk0);

      res_L_e += ai*k*.5*(-M2/k2)*( r00.lo + g2*r00.nlo )*emk0/(1.-emk0);
      res_L_e += ai*k*.5*(-M2/k2)*(rr00.lo +g2*rr00.nlo )*emk0/(1.-emk0);

      res_L_m += ai*k*.5*(-M2/k2)*( r00.lo + g2*r00.nlo )*emk0/(1.-emk0);
      res_L_m += ai*k*.5*(-M2/k2)*(rr00.lo +g2*rr00.nlo )*emk0/(1.-emk0);
    }
    prefactor = 2.*M_PI*pow(T,3.)/(3.*pow(M_PI,3.)*M);
    res_T_e *= prefactor*B_e;
    res_T_m *= prefactor*B_mu;
    res_L_e *= prefactor*B_e;
    res_L_m *= prefactor*B_mu;

    fout << scientific << T
         <<     "    " << muB
         <<     "    " << M*T
         <<     "    " << res_T_e
         <<     "    " << res_T_m
         <<     "    " << res_L_e
         <<     "    " << res_L_m
         << endl;

    //cout << "temp/GeV = " << T << ", M/GeV=" << M << ", k/GeV=" << k << endl;
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();//*/

  return 0;
}


int hydro_table_unintegrated(string mesh_name) {
  int N_kT;
  double res, s, kT_min, kT_max;
  double y = 0.;
  double alpha, B_e, B_mu;

  string fname = "mesh_GJ_nlo.dat";
  fout.open(fname);
  fout << "# Columns: T/GeV, M/GeV, k/GeV, rate" << endl;

  Rho_V rV;
  Rho_00 r00;

  double T,M;
  fin.open("mesh2_alpha.dat");
  //for (int i=0;i<32000;i++) 
  fin.ignore(64,'\n');

  while (!fin.eof()) {
    fin >> T >> M >> k >> alpha;
    fin.ignore(64,'\n');

    cout << "T=" << T << ", M=" << M << ", k=" << k << ", alpha=" << alpha << endl;
    B_e  = B_(m_e*m_e/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.);
    B_mu = B_(m_m*m_m/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.);
    M = M/T;
    k = k/T;
    k0= sqrt( M*M+k*k );

    rV();
    //r00();

    fout << scientific << T
         <<     "    " << M*T
         <<     "    " << k*T
         <<     "    " << B_e*( rV.lo + alpha*4.*M_PI*rV.nlo )*exp(-k0)/(3.*pow(M_PI,3.)*M*M*(1.-exp(-k0)))
         <<     "    " << B_mu*( rV.lo + alpha*4.*M_PI*rV.nlo )*exp(-k0)/(3.*pow(M_PI,3.)*M*M*(1.-exp(-k0)))
         << endl;

    //cout << "temp/GeV = " << T << ", M/GeV=" << M << ", k/GeV=" << k << endl;
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();//*/

  return 0;
}

