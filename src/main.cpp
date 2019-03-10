/*
 *    thermal self-energies, NLO master functions
 *    @author  GJ
 *    @version 0.0
 *
 */
#include "core.hh"
#include <iostream>
double k0, k;
int s[] = {1,1,1};

using namespace std;


int main() {
  master *Rho;
  Rho = _10120(0,0,s);
  cout << "m=" << Rho->m << endl;
  cout << "n=" << Rho->m << endl;
  cout << "s0=" << Rho->s[0] << 
      " ,  s1=" << Rho->s[1] <<
      " ,  s2=" << Rho->s[2] <<
      " ,  s3=" << Rho->s[3] <<
      " ,  s4=" << Rho->s[4] <<
      " ,  s5=" << Rho->s[5] << endl << endl;//*/

  double ope = 48./OOFP;
  cout << (*Rho)(60.99,1.)*ope << endl;
  cout << (*Rho)(100.01,1.)*ope << endl;
  //cout << (*Rho)(2.,1.) << endl;
  //cout << (*Rho)(2.,2.) << endl;
  //cout << k0 << ", " << k << endl;
  //cout << (*Rho)(-1.,3.) << endl;
}
