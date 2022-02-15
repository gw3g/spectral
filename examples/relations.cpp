#include "core.hh"
#include <fstream>
#include <string>

double k0,k;

using namespace std;

int check_1(int,int); 
int check_2(int,int); 

int main(int argc, char *argv[]) {

  if (argc>1) { k0=atof(argv[1]); k=atof(argv[2]); }

  check_2(+1,+1);
  check_2(+1,-1);
  check_2(-1,+1);
  check_2(-1,-1);

}

// RULE 1
int check_1(int sA, int sB)
{
  int S[3];

  S[0] = sA; S[1] = sA; // <-- must be equal!
  S[2] = sB;

  Master *rho_1, *rho_2, *rho_3;
  double      x,      y,      z;

  rho_1 =  _11110(0,0,S);
  rho_2 =  _11110(1,0,S);
  rho_3 =  _11110(0,1,S);

  x = (*rho_1)(k0,k);
  y = (*rho_2)(k0,k);
  z = (*rho_3)(k0,k);
  cout << k0*x-y-2.*z << endl;

  return 0;
}

// RULE 2
int check_2(int sA, int sB)
{
  int S[3];

  S[0] = +1; 
  S[1] = sA;
  S[2] = sB;

  Master *rho_1, *rho_2;
  double      x,      y;

  rho_1 =  _11110(0,0,S);
  rho_2 =  _11110(1,0,S);

  k0 = 5.; k = 1.;
  x = (*rho_1)(k0,k);
  y = (*rho_2)(k0,k);
  cout << k0*x-2.*y << endl;

  return 0;
}
