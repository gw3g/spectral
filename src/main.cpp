#include "core.hh"
#include <iostream>
double k0, k;

using namespace std;


int main() {
  master *Rho;
  int s[] = {1,1,-1};
  Rho = new rho11100(0,2,s);
  cout << "m=" << Rho->m << endl;
  cout << "n=" << Rho->m << endl;
  cout << "s0=" << Rho->s[0] << 
      " ,  s1=" << Rho->s[1] <<
      " ,  s2=" << Rho->s[2] <<
      " ,  s3=" << Rho->s[3] <<
      " ,  s4=" << Rho->s[4] <<
      " ,  s5=" << Rho->s[5] << endl << endl;

  cout << (*Rho)(2.,0.) << endl;
  cout << (*Rho)(2.,1.) << endl;
  cout << (*Rho)(2.,2.) << endl;
  cout << k0 << ", " << k << endl;
  cout << (*Rho)(-1.,3.) << endl;

  cout << "\nThermal\n" << f(.0,-1) << endl;
}
