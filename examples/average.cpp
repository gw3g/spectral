#include "core.hh"
#include "timer.hh"
#include <fstream>
#include <string>

double k0,k;
int S[] = {+1,-1,-1}; // statistics

using namespace std;

ofstream fout;
int elapsed; float percentage;
int Print_k2av();
int V_k2av();

int main() {
  //Print_k2av();
  V_k2av();
}

/*--------------------------------------------------------------------*/

#include <gsl/gsl_sf_bessel.h>
double k2av(double M) {
  // worth noting: k0=3*k for k0~26.
  return 3.*M*gsl_sf_bessel_Kn(3,M)/gsl_sf_bessel_Kn(2,M);
}

/*--------------------------------------------------------------------*/

int V_k2av() {
  int N_M;
  double res, s, M, M_min, M_max;

  Master
    *rho_1;
  double
    _1;

  S[0]=+1; S[1]=+1; S[2]=-1;
  rho_1 =  _11111(2,0,S);


  cout << ":: Creating table ..." << endl << endl;

  string fname = "diag.6.(++-).k2av.dat";
  fout.open(fname);

  fout << "# Columns: M/T, rho*k0/T2, ope (LO), ope (NLO)" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_M=200; 

  M_min=1e-1;
  M_max=1e+2;
  // don't change anything after that.

  s=pow(M_max/M_min,1./(N_M-1));
  M=M_min;

  for (int i=0;i<N_M;i++) { 
    percentage=(float)i/((float)N_M);

    k=sqrt( fabs(k2av(M)) );
    k0=sqrt(M*M+k*k);

    _1 = (*rho_1 )(k0,k);

    fout << scientific  <<  M      // M/T
         << "    "      <<  _1*K2
         << "    "      <<  (rho_1->OPE).T2*K2
         << "    "      <<  (rho_1->OPE).T4*K2
         << endl;

    M*=s;
  }
  cout << endl;
  cout << endl;
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

int Print_k2av() {
  int N_M;
  double res, s, M, M_min, M_max;

  Master
    *rho_b, *rho_bb, *rho_d, *rho_db, *rho_g, *rho_hp, *rho_h_0, *rho_h_1, *rho_j_0, *rho_j_2;
  double
    _b, _bb, _d, _db, _g, _hp, _h_0, _h_1, _j_0, _j_2;

  rho_b   =  _11020(0,0,S);
  rho_bb  =  _10120(0,0,S);
  rho_d   =  _11010(0,0,S);
  rho_db  =  _10110(0,0,S);
  rho_g   =  _11011(0,0,S);
  rho_hp  =  _Star( 0,0,S);
  rho_h_0 =  _11110(0,0,S);
  rho_h_1 =  _11110(0,1,S);
  rho_j_0 =  _11111(0,0,S);
  rho_j_2 =  _11111(2,0,S);


  cout << ":: Creating table ..." << endl << endl;
  
  string fname = "list_k2av.dat";
  fout.open(fname);
  
  fout << "# Columns: M, b, bb, d, db, g, hp, h0, h1, j0, j2" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_M=200; 

  M_min=1e-1;
  M_max=1e+2;
  // don't change anything after that.

  s=pow(M_max/M_min,1./(N_M-1));
  M=M_min;

  for (int i=0;i<N_M;i++) { 
    percentage=(float)i/((float)N_M);

    k=sqrt( fabs(k2av(M)) );
    k0=sqrt(M*M+k*k);

    _b   = (*rho_b )(k0,k)*K2;
    _bb  = (*rho_bb)(k0,k)*K2;
    _d   = (*rho_d )(k0,k);
    _db  = (*rho_db)(k0,k);
    _g   = (*rho_g )(k0,k)*K2;
    _hp  = (*rho_hp)(k0,k);
    _h_0 = (*rho_h_0 )(k0,k)*K2;
    _h_1 = (*rho_h_1 )(k0,k)*k0;
    _j_0 = (*rho_j_0 )(k0,k)*SQR(K2);
    _j_2 = (*rho_j_2 )(k0,k)*K2;



    fout << scientific  <<  M      // M/T
         << "    "      <<  _b
         << "    "      <<  _bb
         << "    "      <<  _d
         << "    "      <<  _db
         << "    "      <<  _g
         << "    "      <<  _hp
         << "    "      <<  _h_0
         << "    "      <<  _h_1
         << "    "      <<  _j_0
         << "    "      <<  _j_2
         << endl;

    M*=s;
  }
  cout << endl;
  cout << endl;
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}
