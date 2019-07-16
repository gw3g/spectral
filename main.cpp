/*
 *    thermal self-energies, NLO Master functions
 *    @author  GJ
 *    @version 0.0
 *
 */
#include "core.hh"
#include <fstream>
#include <string>

double k0, k;
int s[] = {+1,-1,-1};
char s_name[] = {'-','+'};

using namespace std;

ofstream fout;
int Print_k0(Master*,double); int elapsed; float percentage;
void config(Master *rho);

int main() {
  Master *rho;
  rho = _11111(1,1,s);
  config(rho);

  Print_k0(rho,.1);
  Print_k0(rho,1.);
  Print_k0(rho,10.);
  /*double kL = k*(1.-1e-4),
         kR = k*(1.+1e-4);
  double disc= (*rho)(kL,k)*(kL*kL-k*k)
              -(*rho)(kR,k)*(kR*kR-k*k) ;
  cout << " k =  " << k << endl;
  cout << " disc = " << disc << endl;
  cout << " %*64PI*K2 = " << disc*64.*M_PI << endl;
  cout << " %*64*(12/PI)*K2 = " << disc*64.*12/M_PI << endl;//*/
}

void config(Master *rho) {
  // a small function: prints m,n, etc
  cout << "\n:: Configuration of master, TYPE-" << rho->type << endl << endl;
  cout << "m = " << rho->m << endl;
  cout << "n = " << rho->n << endl;
  cout << "s0 = " << rho->s[0] << 
      " ,  s1 = " << rho->s[1] <<
      " ,  s2 = " << rho->s[2] <<
      " ,  s3 = " << rho->s[3] <<
      " ,  s4 = " << rho->s[4] <<
      " ,  s5 = " << rho->s[5] << endl << endl;//*/
}

int Print_k0(Master *rho, double k_curr) {
  int N_k0;
  //percentage;
  double res, s, k0_min, k0_max;
  k=k_curr;

  // filename
  char k_name[20];
  sprintf(k_name,"{k=%.2f}",k);
  string fname = "out/data/diag."
               + to_string(rho->type)
               + string(k_name)
               + ".("+s_name[(rho->s[0]+1)/2]+
                      s_name[(rho->s[1]+1)/2]+
                      s_name[(rho->s[2]+1)/2]+")."
               + to_string(rho->m)+to_string(rho->n)
               + ".dat";

  cout << ":: Creating table for k = " << k <<  " ..." << endl << endl;
  fout.open(fname);
  fout << "# Columns: k0/T, rho, ope (LO), ope (NLO)" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_k0=100; 

  k0_min=1e-2;
  k0_max=1e+2;
  //k0_min = .9*k;
  //k0_max = 1.1*k;
  // don't change anything after that.
  //int o=0; // flag for OPE

  s=pow(k0_max/k0_min,1./(N_k0-1));
  k0=k0_min;

  for (int i=0;i<N_k0;i++) { 
    percentage=(float)i/((float)N_k0);
    /*if (o==1) { 
      res = (*rho)(k0,k);
      fout << k0 << "    " << (rho->OPE)() << endl; 
      k0*=s; continue; }*/
    res = (*rho)(k0,k);
    //if (o==0) { if ( (k0>1.5*k)&&(fabs(1.-res/(rho->OPE)())<1e-2) ) { o=1; }; }
    fout << scientific << k0 << 
                  "    " << res << 
                  "    " << (rho->OPE).T2  <<
                  "    " << (rho->OPE).T4 
                  << endl;
    //fout << k0 << "    " << res << endl;
    k0*=s; 
  }
  cout << endl << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}
