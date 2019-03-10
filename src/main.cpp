/*
 *    thermal self-energies, NLO master functions
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
int s[] = {-1,-1,+1};

using namespace std;

ofstream fout;
int Print_k0(master*,double);
void config(master *Rho);

int main() {
  master *Rho;
  Rho = _11011(0,0,s);
  cout << (*Rho)(1.1,1.) << endl;
  cout << (*Rho)(.9,1.) << endl;
  //cout << k0 << ", " << k << endl;
  //cout << (*Rho)(-1.,3.) << endl;
}

void config(master *Rho) {
  // a small function: prints m,n,etc
  cout << "m=" << Rho->m << endl;
  cout << "n=" << Rho->m << endl;
  cout << "s0=" << Rho->s[0] << 
      " ,  s1=" << Rho->s[1] <<
      " ,  s2=" << Rho->s[2] <<
      " ,  s3=" << Rho->s[3] <<
      " ,  s4=" << Rho->s[4] <<
      " ,  s5=" << Rho->s[5] << endl << endl;//*/
}

int Print_k0(master *Rho, double k_curr) {
  int N_k0;
  double res, s, k0_min, k0_max;
  char *prefix=(char*)"out/data/diagram2_--+_";
  char  suffix[20];
  char  filename[50];

  k=k_curr;
  //signal( SIGALRM, sigalrm_handler );
  //elapsed=0; alarm(10);
  cout << "Creating table for k = " << k <<  " ..." << endl;

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{k=%g}.dat",k);
  strcat(filename,suffix);
  fout.open(filename);

  // Here are some parameters that can be changed:
  N_k0=1000; 
  k0_min=1e-2;
  k0_max=1e2;
  // don't change anything after that.

  s=pow(k0_max/k0_min,1./(N_k0-1));
  k0=k0_min;

  for (int i=0;i<N_k0;i++) { 
    //percentage=(float)i/((float)N_k0);
    fout << k0 << "    " << (*Rho)(k0,k) << endl;
    k0*=s; 
  }
  cout << "Saved to file [" << filename << "]" << endl;
  fout.close();

  return 0;
}
