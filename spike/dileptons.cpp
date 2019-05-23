#include "core.hh"
#include <fstream>
#include <string>

double k0,k;

using namespace std;

ofstream fout;
int print_D(double);
int print_k2av();

int main() {
  print_D(.3);
}

/*--------------------------------------------------------------------*/

#include <signal.h>
#include <unistd.h>
#include <iomanip>
int elapsed;
float percentage=0.;
static void sigalrm_handler( int sig ) {
  elapsed+=1; // monitor progress every 10 sec
  cout.flush();
  cout << "  "<< setw(2) << setfill('0') << elapsed/60 ;
  cout << ":" ;
  cout << setw(2) << setfill('0') << elapsed%60 <<" "; //<< right << 
        //setw(10) << "K = ("<< k0 << "," << k << ")" <<endl;
  cout << '[';
  for (int i=0;i<50;i++) {
    if ((double)i<percentage*50.) {cout << '#';}
    else {cout << '-';};
  }
  cout << "] " << setw(2) << setfill(' ') << (int)(percentage*100.);
  cout << "%" << '\r';
  alarm(1);
}

#include <gsl/gsl_sf_bessel.h>
double k2av(double M) {
  // worth noting: k0=3*k for k0~26.
  return 3.*M*gsl_sf_bessel_Kn(3,M)/gsl_sf_bessel_Kn(2,M);
}

/*--------------------------------------------------------------------*/

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

  Rho_V() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    rho_b =  _10120(0,0,S); // 'vector' channel
    rho_bb=  _11020(0,0,S);
    rho_d =  _10110(0,0,S);
    rho_db=  _11010(0,0,S);
    rho_g =  _11011(0,0,S);
    rho_h =  _11110(0,0,S);
    rho_hp=  _Star( 0,0,S);
    rho_j =  _11111(0,0,S);

  };
  void operator ()() {

    //lo = 2.*Nc*K2*lga( cosh(.5*kp)/cosh(.5*km) )/k*OOFP;
    lo = 2.*Nc*K2*(.5+lga( (1.+exp(-kp))/(1.+exp(-km)) ))/k*OOFP;
    nlo= lo*3.*cF*SQR(OOFP);

    _b = (*rho_b )(k0,k)*K2;
    _bb= (*rho_bb)(k0,k)*K2;
    _d = (*rho_d )(k0,k);
    _db= (*rho_db)(k0,k);
    _g = (*rho_g )(k0,k)*K2;
    _h = (*rho_h )(k0,k);
    _hp= (*rho_hp)(k0,k);
    _j = (*rho_j )(k0,k);

    nlo +=
    8.*Nc*cF*( 2.*(_b-_bb+_d-_db) + 1.5*_g - 2.*(_h+_hp) + _j );

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
    _b_0, _bb_0,
    _b_1, _bb_1,
    _b_2, _bb_2,
    _g,
    _h_0, _h_1,
    _hp,
    _j_0, _j_2;

  Rho_00() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    rho_b_0 =  _10120(0,0,S); // notation: rho_<tag>_<power of ..>
    rho_bb_0=  _11020(0,0,S);
    rho_b_1 =  _10120(1,0,S);
    rho_bb_1=  _11020(1,0,S);
    rho_b_2 =  _10120(2,0,S);
    rho_bb_2=  _11020(2,0,S);
    rho_g   =  _11011(0,0,S);
    rho_h_0 =  _11110(0,0,S);
    rho_h_1 =  _11110(0,1,S);
    rho_hp  =  _Star( 0,0,S);
    rho_j_0 =  _11111(0,0,S);
    rho_j_2 =  _11111(2,0,S);

  };
  void operator ()() {

    //lo = 2.*Nc*K2*lga( cosh(.5*kp)/cosh(.5*km) )/k*OOFP;
    lo = 2.*Nc*K2*(.5+lga( (1.+exp(-kp))/(1.+exp(-km)) ))/k*OOFP;
    nlo= lo*3.*cF*SQR(OOFP);

    _b_0 = (*rho_b_0 )(k0,k)*K2;
    _bb_0= (*rho_bb_0)(k0,k)*K2;
    _b_1 = (*rho_b_1 )(k0,k)*k0;
    _bb_1= (*rho_bb_1)(k0,k)*k0;
    _b_2 = (*rho_b_1 )(k0,k);
    _bb_2= (*rho_bb_1)(k0,k);
    _g   = (*rho_g   )(k0,k)*k*k;
    _h_0 = (*rho_h_0 )(k0,k);   // TODO: check norm of h, j
    _h_1 = (*rho_h_1 )(k0,k)*k0/K2;
    _hp  = (*rho_hp)(k0,k);
    _j_0 = (*rho_j_0 )(k0,k)*(k0*k0+k*k)/K2;
    _j_2 = (*rho_j_2 )(k0,k)*k0*k0/K2;

    nlo +=
    8.*Nc*cF*( 2.*(_b_0-_bb_0-4.*(_b_1-_bb_1)+4.*(_b_2-_bb_2)) 
             - k*k*_g - 2.*(_h_0+_hp) + 8.*_h_1 - _j_0 - 4.*_j_2 );

  };
};


int print_D(double k_curr) {
  int N_k0;
  double res, s, k0_min, k0_max;
  k=k_curr;

  // filename
  char k_name[20];
  sprintf(k_name,"{k=%.2f}",k);
  string fname = "spike/NLO_rhoV_"
               + string(k_name)
               + ".dat";

  cout << ":: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V _R;
  fout.open(fname);
  fout << "# Columns: k0/T, -rhoV_LO/T2, -rhoV_NLO/(g2*T2)" << endl;
  fout << "# ( k=" << k << " )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_k0=600; 

  k0_min=1e-2;
  k0_max=1e+2;
  // don't change anything after that.

  s=pow(k0_max/k0_min,1./(N_k0-1));
  k0=k0_min;

  for (int i=0;i<N_k0;i++) { 
  //while (k0<10.) {
    percentage=(float)i/((float)N_k0);
    _R();
    //cout << k0 << "    " << res << endl;
    fout << k0     << "    "
         << _R.lo  << "    "
         << _R.nlo << endl;

    k0*=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

int print_k2av() {
  int N_M;
  double res, s, M, M_min, M_max;

  cout << ":: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V _R;
  string fname = "spike/list.Mr_k2av.dat";
  fout.open(fname);
  fout << "# Columns: b, bb, d, db, g, h, hp, j" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_M=400; 

  M_min=1e-1;
  M_max=1e+2;
  // don't change anything after that.

  s=pow(M_max/M_min,1./(N_M-1));
  M=M_min;

  for (int i=0;i<N_M;i++) { 
    percentage=(float)i/((float)N_M);

    k=sqrt( fabs(k2av(M)) );
    k0=sqrt(M*M+k*k);
    _R();

    fout << setprecision(9) << M 
      << "    " <<  _R._b
      << "    " <<  _R._bb
      << "    " <<  _R._d
      << "    " <<  _R._db
      << "    " <<  _R._g
      << "    " <<  _R._h
      << "    " <<  _R._hp
      << "    " <<  _R._j
     << endl;
    M*=s;
  }
  cout << endl;
  cout << endl;
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}
