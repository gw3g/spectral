#include <cmath>
#include <iostream>
using namespace std;

#define sgn(x) (double) ((x>0)-(x<0))
#define lga(x) log(fabs(x))
#define clp(x,a,b) max(a,min(b,x))
#define OOFP 0.0795774715459476678844418816863

template<typename T>
inline T SQR(const T x) {return x*x;}

template<typename T>
inline T CUBE(const T x) {return x*SQR(x);}

// thermal distribution functions
double  f(double,int),
       fb(double,int);

double I(int,int); // tadpole
double ff(double , double, int , int);
double fff(double , double, double, int , int , int *);

void print_integrand(int m, int n, int s[3]); // for checking

/*--------------------------------------------------------------------*/

#define kp ((k0+k)*.5)
#define km ((k0-k)*.5)
#define K2 (k0*k0-k*k)
extern double k0, k; // external 4-momentum

struct Master {
  int s[6]; // statistical configuration
  int m, n; // p^n.q^m
  virtual double eval()=0;
  Master(int _m, int _n, int _s[3])
    : m(_m), n(_n) {
        // provided:
        s[0] = _s[0]; s[1] = _s[1]; s[2] = _s[2];
        // derived:
        s[3] = _s[0]*_s[1]*_s[2];
        s[4] = _s[0]*_s[1];
        s[5] = _s[0]*_s[2];
    }
  double operator ()(double _k0, double _k) { 
    k0=_k0; k=_k; // reset K
    return eval();
  };
};

/*--------------------------------------------------------------------*/

// type.I
Master* _01020(int, int, int *);
Master* _00120(int, int, int *);

// type.II
Master* _11010(int, int, int *);
Master* _11020(int, int, int *);
Master* _10120(int, int, int *);

// type.III
Master* _11011(int, int, int *);

// type.IV
Master* _11100(int, int, int *);

// type.V
Master* _11110(int, int, int *);
Master*  _Star(int, int, int *);

// type.VI
Master* _11111(int, int, int *);

