/*
 *    thermal self-energies, NLO Master functions
 *    @author  GJ
 *    @version 0.2
 *
 */
#include "core.hh"
#include "timer.hh"
#include <fstream>
#include <string>

double k0, k;
double mu_q = 0.;
bool  CHEM_POT = false;
double MOT1, MOT2, MOT3, MOT4, MOT5;// mu over T
int s[3];// = {+1,-1,-1};
char s_name[] = {'-','+'};

using namespace std;

ofstream fout;
ifstream fin;
int Print_k0(Master*,double); int elapsed; float percentage;
void config(Master *rho);

int main(int argc, char *argv[]) {
  char cc;
  char *str_in;

  Master *rho;
  fin.open("config");

  char _t[6]; fin >> _t; // abcde
  char _s[4]; fin >> _s; // (±±±)
  int _m, _n; fin >> _m; fin >> _n;

  cout << "\n:: Input [config]\n" << endl; // regurgitate
  cout << _t       << endl;
  cout << _s       << endl;
  cout << _m       << endl;
  cout << _n       << endl;

  s[0] = (_s[0]=='+') ? +1 : -1;
  s[1] = (_s[1]=='+') ? +1 : -1;
  s[2] = (_s[2]=='+') ? +1 : -1;

       if (string(_t)=="01020\0") { rho = _01020(_m,_n,s); }
  else if (string(_t)=="00120\0") { rho = _00120(_m,_n,s); }
  else if (string(_t)=="11020\0") { rho = _11020(_m,_n,s); }
  else if (string(_t)=="10120\0") { rho = _10120(_m,_n,s); }
  else if (string(_t)=="11010\0") { rho = _11010(_m,_n,s); }
  else if (string(_t)=="10110\0") { rho = _10110(_m,_n,s); }
  else if (string(_t)=="11011\0") { rho = _11011(_m,_n,s); }
  else if (string(_t)=="11100\0") { rho = _11100(_m,_n,s); }
  else if (string(_t)=="11110\0") { rho = _11110(_m,_n,s); }
  else if (string(_t)=="Star\0")  { rho = _Star(0,0,s); }
  else if (string(_t)=="11111\0") { rho = _11111(_m,_n,s); }
  else    { cout << " Invalid! " << endl; return -1; }

  --argc;
  while ( (++argv)[0] && argv[0][0]=='-' ) { // parse arguments
    --argc;
    while ((cc = *(++argv[0]))) {
      switch (cc) {
        case 'k': // 3-momentum
          if (argc==0) { k=1.0; }
          else { 
            str_in=(*++argv);
            k = atof(str_in);
          }
          Print_k0(rho,k);
          break;
        case 'm':
          CHEM_POT = true;
          if (argc==0) { mu_q=.0; }
          else { 
            str_in=(*++argv);
            mu_q = atof(str_in);
          }
          MOT1 = mu_q; MOT2 = -mu_q;
          MOT3 = - MOT1 - MOT2; MOT4 = -MOT1; MOT5 = -MOT2;
          Print_k0(rho,1.);
          break;
        case 'c':
          config(rho);
          break;
        case 'h':
          cerr << "\n:: Usage, './bin/rho' ..."                         << endl << endl;
          cerr << "         [-k momentum] (in units of T, default=1.0)" << endl;
          cerr << "         [-c] (show master config)"                  << endl;
          cerr << "         [-h] (prints this help)"                    << endl;
          break;
      }
    }
  }
  cout << endl; return 0;
}

void config(Master *rho) {
  // a small function: prints m,n, etc
  cout << "\n:: Master, TYPE-" << rho->type;
  if (rho->type == 7) { cout << " (Star)"; }
  cout << endl << endl;
  cout << "m = " << rho->m << endl;
  cout << "n = " << rho->n << endl;
  cout << "s0 = " << rho->s[0] << 
      " ,  s1 = " << rho->s[1] <<
      " ,  s2 = " << rho->s[2] <<
      " ,  s3 = " << rho->s[3] <<
      " ,  s4 = " << rho->s[4] <<
      " ,  s5 = " << rho->s[5] << endl;//*/
}

int Print_k0(Master *rho, double k_curr) {
  int N_k0;
  //percentage;
  double res, s, k0_min, k0_max;
  k=k_curr;

  // filename
  char k_name[20];
  char mu_name[20];
  sprintf(k_name,"{k=%.2f}",k);
  sprintf(mu_name,"{mu=%.2f}",mu_q);
  string fname = "out/data/diag."
               + to_string(rho->type)
               + string(k_name)
               + ".("+s_name[(rho->s[0]+1)/2]+
                      s_name[(rho->s[1]+1)/2]+
                      s_name[(rho->s[2]+1)/2]+")."
               + to_string(rho->m)+to_string(rho->n);
  if (CHEM_POT) {fname = fname+"."+string(mu_name); }
  fname = fname + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  fout.open(fname);
  fout << "# Columns: k0/T, rho, ope (LO), ope (NLO)" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_k0=20; 

  k0_min=1e+1;
  k0_max=1e+2;

  // don't change anything after that.
  // int o=0; // flag for OPE

  s=pow(k0_max/k0_min,1./(N_k0-1));
  k0=k0_min;

  for (int i=0;i<N_k0;i++) { 
    percentage=(float)i/((float)N_k0);
    res = (*rho)(k0,k);

    fout << scientific << k0 
         <<     "    " << res
         <<     "    " << (rho->OPE).T2
         <<     "    " << (rho->OPE).T4 
         << endl;

    k0*=s; 
  }

  cout << "  "<< setw(2) << setfill('0') << elapsed/60 ;
  cout << ":" ;
  cout << setw(2) << setfill('0') << elapsed%60 <<" ";
  cout << '[';
  for (int i=0;i<50;i++) { cout << "#"; }
  cout << "] " << "100%";

  cout << endl << endl << ":: Saved to file [" << fname << "]\n";
  fout.close();

  return 0;
}
