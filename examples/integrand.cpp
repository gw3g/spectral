
#include "core.hh"
#include <fstream>


double k0,k;

using namespace std;
int elapsed; float percentage;

void print_integrand(int m, int n, int s[3]) {
  ofstream fout;

  Master *rho;
  rho =  _11110(0,0,s);

  fout.open("test_integrand.dat");

  double x=-.0125*.5;
  double y;
  for (int i=0;i<80;i++) {
    x+=.025*.5;
    y=-.0125*.5;
    for (int j=0;j<80;j++) {
      y+=.025*.5;
      fout  << x << "    " 
            << y << "    " 
            << (rho->integrand)(x,y)  
            << endl; 
    }
  }
  fout.close(); 
}//*/

int main(int argc, char *argv[]) {

  int S[3] = {+1,+1,+1};
  if (argc>1) { k0=atof(argv[1]); k=atof(argv[2]); }
  print_integrand(0,0,S);

}
