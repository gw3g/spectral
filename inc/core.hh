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

double psi0(int sA, int sB); // moments
double psi1(int sA, int sB);
double psi2(int sA, int sB);

void print_integrand(int m, int n, int s[3]); // for checking

/*--------------------------------------------------------------------*/

#define kp ((k0+k)*.5)
#define km ((k0-k)*.5)
#define K2 (k0*k0-k*k)
extern double k0, k; // external 4-momentum

struct Expand {
  double T0, T2, T4; // coeffs of T^n
  double operator()() { return T2+T4; };
};

struct Master {
  int s[6]; // statistical configuration
  int m, n; // p^n.q^m
  int type; // = {1..6}
  Expand OPE;
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
    // dimension = 2(a+b+c+d+e)-m-n-6
  };
};

/*--------------------------------------------------------------------*/

// type.I
Master* _01020(int, int, int *);
Master* _00120(int, int, int *);

// type.II
Master* _11010(int, int, int *);
Master* _10110(int, int, int *);
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

/*--------------------------------------------------------------------*/
// for output while running:
#include <signal.h>
#include <unistd.h>
#include <iomanip>
extern int elapsed;
extern float percentage;

static void sigalrm_handler( int sig ) {
  elapsed+=1; // monitor progress every 10 sec
  cout.flush();
  cout << "  "<< setw(2) << setfill('0') << elapsed/60 ;
  cout << ":" ;
  cout << setw(2) << setfill('0') << elapsed%60 <<" "; //<< right << 
  //cout << setw(10) << "K = ("<< setprecision(2) << k0 << "," << k << ")" << "  ";
  cout << '[';
  for (int i=0;i<50;i++) {
    if ((double)i<percentage*50.) {cout << '#';}
    else {cout << '-';};
  }
  cout << "] " << setw(2) << setfill(' ') << (int)(percentage*100.);
  cout << "%" << '\r';
  alarm(1);
}

/*--------------------------------------------------------------------*/
