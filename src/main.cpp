/*
 *    thermal self-energies, NLO Master functions
 *    @author  GJ
 *    @version 0.0
 *
 */
#include "core.hh"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
double k0, k;
int s[] = {1,1,1};

using namespace std;

ofstream fout;
int Print_k0(Master*,double);
void config(Master *rho);
char PM(int i) { if (i<0) return '-'; else return '+'; }

int main() {
  Master *rho;
  rho = _11111(0,0,s);

  //k0 = 60.5; k = 1.0;
  //print_integrand(0,0,s);
  //Print_k0(rho,.004);
  //Print_k0(rho,.1);
  Print_k0(rho,1.);
  //Print_k0(rho,10.);
  //cout << k0 << ", " << k << endl;
}

void config(Master *rho) {
  // a small function: prints m,n,etc
  cout << "m =" << rho->m << endl;
  cout << "n =" << rho->m << endl;
  cout << "s0=" << rho->s[0] << 
      " ,  s1=" << rho->s[1] <<
      " ,  s2=" << rho->s[2] <<
      " ,  s3=" << rho->s[3] <<
      " ,  s4=" << rho->s[4] <<
      " ,  s5=" << rho->s[5] << endl << endl;//*/
}

int Print_k0(Master *rho, double k_curr) {
  int N_k0;
  double res, s, k0_min, k0_max;
  char *prefix=(char*)"out/data/diagram6_+++_00";
  char  suffix[20];
  char  filename[50];

  k=k_curr;
  //signal( SIGALRM, sigalrm_handler );
  //elapsed=0; alarm(10);
  cout << "Creating table for k = " << k <<  " ..." << endl;
  //string si = (string)PM(rho->s[0])+PM(rho->s[1])+PM(rho->s[2]);

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{k=%g}.dat",k);
  strcat(filename,suffix);
  fout.open(filename);

  // Here are some parameters that can be changed:
  N_k0=50; 
  k0_min=5e1;
  k0_max=1e2;
  // don't change anything after that.

  s=pow(k0_max/k0_min,1./(N_k0-1));
  k0=k0_min;

  for (int i=0;i<N_k0;i++) { 
    //percentage=(float)i/((float)N_k0);
    res = (*rho)(k0,k);
    cout << k0 << "    " << res << endl;
    fout << k0 << "    " << res << endl;
    k0*=s; 
  }
  cout << "Saved to file [" << filename << "]" << endl;
  fout.close();

  return 0;
}
