#include "core.hh"
#include "quad.hh"
#include "map.hh"
#include <fstream>
#include <string>

double k0,k;
double mu_q = 0.;
bool  CHEM_POT = false;
double MOT1, MOT2, MOT3, MOT4, MOT5;// mu over T
double tol = 1e0; // absolute tolerance

using namespace std;

ofstream fout;
ifstream fin;
int Print_D(double); int elapsed; float percentage;
int Print_D(double,double); 
void D(double,double);
int hydro_table();
int hydro_table2();
int hydro_table_T_L();

int main(int argc, char *argv[]) {

  // Enter: k0, k
  if (argc>1) { k0=atof(argv[1]); k=atof(argv[2]); }

  // Return: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
  //D(k0,k);

  /* for lattice comparisons: */
  //nf=0
  Print_D(1.*2.*M_PI/3.,0.);     // T=1.1Tc
  Print_D(1.*2.*M_PI/3.,1.);     // T=1.1Tc
  Print_D(1.*2.*M_PI/3.,2.);     // T=1.1Tc
  //Print_D(3.*7.*M_PI/12.);    // T=1.3Tc
  // nf=2
  //Print_D(sqrt(1.)*M_PI/2.); // T=1.2Tc
  //Print_D(1.5*M_PI);         // T=1.2Tc
  //Print_D(.5*M_PI*sqrt(14)); // T=1.2Tc
  //
  /* for hydro run: */
  //hydro_table_T_L();

  /* Others: */
  //Print_D(.02);
  //Print_D(1.);
  //Print_D(1.5);

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

  // notation à la 1310.0164
  Master
    *rho_b, *rho_bb, *rho_d, *rho_db, *rho_g, *rho_h, *rho_hp, *rho_j;
  double
    _b, _bb, _d, _db, _g, _h, _hp, _j;

  size_t limit = 5e2;
  quad *wsp1;
  quad *wsp2;

  Rho_V() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    rho_b =  _11020(0,0,S); // 'vector' channel
    rho_bb=  _10120(0,0,S);
    rho_d =  _11010(0,0,S);
    rho_db=  _10110(0,0,S);
    rho_g =  _11011(0,0,S);
    rho_h =  _11110(0,0,S);
    rho_hp=  _Star( 0,0,S);
    rho_j =  _11111(0,0,S);

    wsp1 = new quad(limit); // prepare workspace quadrature
    wsp2 = new quad(limit);
    gsl_set_error_handler_off();

  };
  void operator ()() {

    lo = +Nc*K2*psi0(-1,-1,MOT1)*OOFP; // note sign convention
    //lo = Nc*K2*OOFP; // large-K2
    nlo= lo*3.*cF*SQR(OOFP);

    // Quadrature step! --
    double res, err;

    auto outer = make_gsl_function( [&](double x) 
    {
      double inner_result, inner_abserr;
      auto inner = make_gsl_function( [&](double y) {
            return (  2.*(rho_h ->integrand)(x,y)
                     +2.*(rho_hp->integrand)(x,y)*(kp/km)
                     -   (rho_j ->integrand)(x,y)*(kp/km)
                   )/SQR(kp);
          } );
      gsl_integration_qag( inner, .0,1., tol, 0,
                          limit, 6, *wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, .0,1., tol*2.,  0,
                        limit, 6, *wsp2, &res, &err  );//*/

    // Simpler masters --
    _b = (*rho_b )(k0,k)*K2;
    _bb= (*rho_bb)(k0,k)*K2;
    _d = (*rho_d )(k0,k);
    _db= (*rho_db)(k0,k);
    _g = (*rho_g )(k0,k)*K2;
    //_h = (*rho_h )(k0,k)*K2;
    //_hp= (*rho_hp)(k0,k);
    //_j = (*rho_j )(k0,k)*SQR(K2);

    nlo -=
    //8.*Nc*cF*( 2.*(_b-_bb+_d-_db) - 1.5*_g + 2.*(_h+_hp) - _j );
    8.*Nc*cF*( 2.*(_b-_bb+_d-_db) - 1.5*_g + res*CUBE(OOFP)*SQR(kp) );

  };
};

struct Rho_00
{
  double lo, nlo;
  int S[3];
  double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc);

  Master
    *rho_b_0, *rho_bb_0,
    *rho_b_1, *rho_bb_1,
    *rho_b_2, *rho_bb_2,
    *rho_g,
    *rho_h_0, *rho_h_1,
    *rho_hp,
    *rho_j_0, *rho_j_2;
  double
    _b_1, _bb_1,
    _b_2, _bb_2,
    _g,
    _h_0, _h_1,
    _hp,
    _j_0, _j_2;

  size_t limit = 5e2;
  quad *wsp1;
  quad *wsp2;

  Rho_00() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    //rho_b_0 =  _11020(0,0,S); // notation: rho_<tag>_<power of ..>
    //rho_bb_0=  _10120(0,0,S);
    //rho_b_1 =  _11020(1,0,S);
    //rho_bb_1=  _10120(1,0,S);
    //rho_b_2 =  _11020(2,0,S);
    //rho_bb_2=  _10120(2,0,S);
    rho_g   =  _11011(0,0,S);
    rho_h_0 =  _11110(0,0,S);
    rho_h_1 =  _11110(0,1,S);
    rho_hp  =  _Star( 0,0,S);
    rho_j_0 =  _11111(0,0,S);
    rho_j_2 =  _11111(2,0,S);

    wsp1 = new quad(limit);
    wsp2 = new quad(limit);
    gsl_set_error_handler_off();

  };
  void operator ()() {

    lo = -2.*Nc*( k0*k0*(psi1(-1,-1,MOT1)-psi2(-1,-1,MOT1))-.25*K2*psi0(-1,-1,MOT1) )*OOFP;
    nlo = -cF*Nc*( k*k*psi0(-1,-1,MOT1) )*CUBE(OOFP);
    //nlo = -cF*Nc*( k*k )*CUBE(OOFP); // large-K2

    // Quadrature step! --
    double res, err;

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
      gsl_integration_qag( inner, .0,1., tol, 0,
                          limit, 6, *wsp1, &inner_result, &inner_abserr );
      return inner_result;
    } );
    gsl_integration_qag( outer, .0,1., tol*2., 0,
                       limit, 6, *wsp2, &res, &err  );//*/

    // Simpler master(s) --
    //_b_0 = (*rho_b_0 )(k0,k)*K2;
    //_bb_0= (*rho_bb_0)(k0,k)*K2;
    //_b_1 = (*rho_b_1 )(k0,k)*k0;
    //_bb_1= (*rho_bb_1)(k0,k)*k0;
    //_b_2 = (*rho_b_2 )(k0,k);
    //_bb_2= (*rho_bb_2)(k0,k);
    _g   = (*rho_g   )(k0,k)*k*k;
    //_h_0 = (*rho_h_0 )(k0,k)*K2;
    //_h_1 = (*rho_h_1 )(k0,k)*k0;
    //_hp  = (*rho_hp  )(k0,k);
    //_j_0 = (*rho_j_0 )(k0,k)*(k0*k0+k*k)*K2;
    //_j_2 = (*rho_j_2 )(k0,k)*K2;

    nlo -=
    4.*Nc*cF*( //2.*(_b_0-_bb_0-4.*(_b_1-_bb_1)+4.*(_b_2-_bb_2)) // =0 
              //+ _g + 2.*(_h_0+_hp) - 8.*_h_1 + _j_0 - 4.*_j_2 );
             + _g + res*CUBE(OOFP)*SQR(kp) );

  };
};

/*--------------------------------------------------------------------*/
// evaluation

void D(double k0_curr, double k_curr) {
  k0=k0_curr; k=k_curr;

  Rho_V   rV; // assign ptrs
  Rho_00 r00;

  rV(); //r00();


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
  fname = fname + ".2.dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV, rV2;
  Rho_00 r00;
  //fout.open(fname);
  cout << 
  "# Columns: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
       << endl
       << "# ( k=" << k << " )" << endl;

  //signal( SIGALRM, sigalrm_handler );
  //elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_k0=400; 

  k0_min=1e-4;
  k0_max=10.;
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  s=(k0_max-k0_min)/((double)N_k0-1.);
  //s = 1e-1;
  k0+=s;
  k0=k0_min;

  int i=0;
  //for (int i=0;i<N_k0;i++) { 
  while (k0<k0_max) {
    //percentage=(float)i/((float)N_k0);
    percentage = (k0-k0_min)/(k0_max-k0_min);

    MOT1 = mu; MOT2 = -mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV();
    MOT1 = -mu; MOT2 = mu;
    MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
    rV2();
    //r00();

    fout << scientific << k0        // k0/T
         <<     "    " << .5*(rV.lo+rV2.lo)     // leading-order: rho_V ,
         <<     "    " << 0.//r00.lo    //                rho_00
         <<     "    " << .5*(rV.nlo+rV2.nlo)    // next-to-LO   : rho_V ,
         <<     "    " << 0.//r00.nlo   //                rho_00
        << endl;

    k0+=s; 
  }
  //cout << endl << ":: Saved to file [" << fname << "]" << endl;
  //fout.close();

  return 0;
}

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

// TO BE TIDIED UP!!
int hydro_table() { //
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
  fin.open("mesh_alpha2.dat");
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


int hydro_table_T_L() { //
  int N_kT;
  double res_T_e, res_T_m, 
         res_L_e, res_L_m,
         s, kT_min=0., kT_max=15.;
  double y = 0.;
  double alpha, B_e, B_mu;
  double ai, ti;

  string fname = "mesh_GJ_kT_integrated_T_L.dat";
  fout.open(fname);
  fout << "# Columns: T/GeV, M/GeV, rate_T(e+e-), rate_T(mu+mu-), rate_L(e+e-), rate_L(mu+mu-)" << endl;

  Rho_V rV;
  Rho_00 r00;

  double T,M;
  fin.open("mesh_alpha2.dat");
  //for (int i=0;i<32000;i++) 
  fin.ignore(64,'\n');

  while (!fin.eof()) {
    fin >> T >> M >> alpha;
    fin.ignore(64,'\n');

    cout << "T=" << T << ", M=" << M << ", alpha=" << alpha << endl;
    B_e  = B_(m_e*m_e/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.); // units:
    B_mu = B_(m_m*m_m/(M*M))*(2./3.)*pow(hbarc,-4.)*pow(alpha_em,2.); // GeV^-4.fm^-4

    M = M/T;

    res_T_e = 0.;
    res_T_m = 0.;
    res_L_e = 0.;
    res_L_m = 0.;

    for (int i=0;i<32;i++) {
      ai  = .5*G32pt[i][0]; ti = .5*(G32pt[i][1]+1.);
      ai *= kT_max-kT_min;
      k   = kT_min*(1.-ti) + kT_max*ti;// = k/T
      k0  = sqrt( M*M+k*k );

      rV();
      r00();

      res_T_e += ai*k*.5*( rV.lo + r00.lo*(M*M/(k*k)) + alpha*4.*M_PI*( rV.nlo + r00.nlo*(M*M/(k*k)) ) )*exp(-k0)/(1.-exp(-k0));
      res_T_m += ai*k*.5*( rV.lo + r00.lo*(M*M/(k*k)) + alpha*4.*M_PI*( rV.nlo + r00.nlo*(M*M/(k*k)) ) )*exp(-k0)/(1.-exp(-k0));
      res_L_e += ai*k*(-M*M/(k*k))*( r00.lo + alpha*4.*M_PI*r00.nlo )*exp(-k0)/(1.-exp(-k0));
      res_L_m += ai*k*(-M*M/(k*k))*( r00.lo + alpha*4.*M_PI*r00.nlo )*exp(-k0)/(1.-exp(-k0));
    }
    res_T_e *= 2.*M_PI*pow(T,3.)*B_e/(3.*pow(M_PI,3.)*M);
    res_T_m *= 2.*M_PI*pow(T,3.)*B_mu/(3.*pow(M_PI,3.)*M);
    res_L_e *= 2.*M_PI*pow(T,3.)*B_e/(3.*pow(M_PI,3.)*M);
    res_L_m *= 2.*M_PI*pow(T,3.)*B_mu/(3.*pow(M_PI,3.)*M);

    fout << scientific << T
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



int hydro_table2() { //
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
  fin.open("mesh_alpha.dat");
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

