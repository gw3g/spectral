#define sgn(x) (double) ((x>0)-(x<0))
#define lga(x) log(fabs(x))
#define clp(x,a,b) max(a,min(b,x))

// thermal distribution functions
double  f(double,int),
       fb(double,int);

/*--------------------------------------------------------------------*/

#define kp ((k0+k)*.5)
#define km ((k0-k)*.5)
#define K2 (k0*k0-k*k)
extern double k0, k; // external 4-momentum

struct master {
  int s[6]; // statistical configuration
  int m, n; // p^n.q^m
  virtual double eval()=0;
  master(int _m, int _n, int _s[3])
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
//
//
master* _11010(int, int, int *);

