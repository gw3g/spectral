#include "core.hh"
#include <fstream>
#include <string>

using namespace std;

ofstream fout;
int print_D(double);

int main() {
  print_D(1.);
}

struct Dilepton {
  // notation a la 1310.0164
  int S[3];
  double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc);
  Master *rho_b;
  Master *rho_bb;
  Master *rho_d;
  Master *rho_db;
  Master *rho_g;
  Master *rho_h;
  Master *rho_hp;
  Master *rho_j;

  Dilepton() {

    S[0] = +1; S[1] = - 1; S[2] = -1;

    rho_b =  _10120(0,0,S);
    rho_bb=  _11020(0,0,S);
    rho_d =  _10110(0,0,S);
    rho_db=  _11010(0,0,S);
    rho_g =  _11011(0,0,S);
    rho_h =  _11110(0,0,S);
    rho_hp=  _Star( 0,0,S);
    rho_j =  _11111(0,0,S);

  };
  double operator ()() {

    double _b = (*rho_b )(k0,k),
           _bb= (*rho_bb)(k0,k),
           _d = (*rho_d )(k0,k),
           _db= (*rho_db)(k0,k),
           _g = (*rho_g )(k0,k),
           _h = (*rho_h )(k0,k),
           _hp= (*rho_hp)(k0,k),
           _j = (*rho_j )(k0,k);

    return 8.*Nc*cF*(
                      2.*( K2*(_b-_bb)+_d-_db) + 1.5*K2*_g - 2.*(_h+_hp) + _j 
                    );
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

  cout << "Creating table for k = " << k <<  " ..." << endl;
  Dilepton _R;
  fout.open(fname);

  // Here are some parameters that can be changed:
  N_k0=20; 

  k0_min=1e-1;
  k0_max=1e+1;
  s = 1e-2;
  //k0_min = .9*k;
  //k0_max = 1.1*k;
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  k0=k0_min;

  //for (int i=0;i<N_k0;i++) { 
  while (k0<10.) {
    //percentage=(float)i/((float)N_k0);
    res = _R() ;
    cout << k0 << "    " << res << endl;
    fout << k0 << "    " << res << endl;
    k0+=s; 
  }
  cout << "Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}
