#define sgn(x) (double) ((x>0)-(x<0))
#define lga(x) log(fabs(x))
#define clp(x,a,b) max(a,min(b,x))
#define OOFP 0.0795774715459476678844418816863

// thermal distribution functions
double  f(double,int),
       fb(double,int);
double I(int,int); // tadpole

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

// Type.I
master* _01020(int, int, int *);
master* _00120(int, int, int *);

// Type.II
master* _11010(int, int, int *);
master* _11020(int, int, int *);
master* _10120(int, int, int *);

// Type.III
master* _11011(int, int, int *);

