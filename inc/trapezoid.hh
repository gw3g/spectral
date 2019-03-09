#include <math.h>
/*
 * 2nd Euler-Maclaurin sum formula,
 * (extended midpoint rule, OPEN)
 *
 */
template <class F>
struct integrate { // func(x) over [0,1]
  double res;
  int n; // mesh resolution ~ 3^n
  F &funk;
  integrate(F &_func) : funk(_func) {n=0;}
  double next() {
    int i,j;
    double x, temp, del, ddel;
    n++;
    if (n==1) { // starting estimate
      return res=func(.5);
    } else { // triple number of steps to benefit from
             // all the previous function evaluations.
      temp=.0;
      for (i=1,j=1;j<n-1;j++) i*=3; // 3^{n-2}
      del=1./(3.*(double)i);
      x=.5*del;
      ddel=del+del; // new points alternate in spacing...
      for (j=0;j<i;j++) {
        temp += func(x);
        x += ddel;  // ie, skip x+=del to not double count
        temp += func(x);
        x += del;
      }
      res = ( res+temp/(double)i )/3.; // refined estimate
      return res;
    }
  }
  // identity mapping, placeholder for derived rules
  virtual double func(double x) {return funk(x);}
};

/*--------------------------------------------------------------------*/
// implement variable changes:
//
//  map_Linear,   [0,1] -> [a,b]
//  map_Finite,   [0,1] -> [a,b] for endpoint singularities
//  map_SemiInf,  [0,1] -> [a,+inf] (or [-inf,a] if a<0)
//
//
struct map {
  double a, b, tmax=4.;
  virtual double operator ()(double t) = 0; // t is the new var
};

template <typename F>
struct Finite : map {
  //double a, b, tmax=4.;
  F &funk;
  double operator ()(double t) {
    double dxdt, del, Q = exp(-2.*sinh(t*tmax));
    dxdt = tmax*2.*(b-a)*Q*cosh(t*tmax)/( (1.+Q)*(1.+Q) );
    del = (b-a)*Q/(1.+Q);
    return ( // trick to avoid cancellation errors
        funk( a+del , b-a-del )+
        funk( b-del , del ) 
        )*dxdt;
  }
  Finite(F &_func, double _a, double _b) :
    funk(_func) { a=_a; b=_b; }
};

template <typename F>
struct SemiInf : map {
  //double a, tmax=4.;
  F &funk;
  double operator ()(double t) {
    double x, dxdt, Q = exp(2.*sinh(((1.-2.*t)*tmax)));
    x = a*(1.+Q);
    dxdt = tmax*4.*a*Q*cosh(((1.-2.*t)*tmax));
    return funk(x,0)*dxdt;
  }
  SemiInf(F &_func, double _a) :
    funk(_func) { a=_a; }
};

/*--------------------------------------------------------------------*/
// runner, eps is relative
//#include <iostream>
//using namespace std;

template <typename F>
double go(integrate<F> &I, double eps=1.e-6) {
  const int iMax=10; // max number of iterations
  double res, old;
  for (int i=0;i<iMax;i++) {
    res = I.next(); //cout << res << endl;
    if (i>5)
      if (abs(res-old)<eps*abs(old)||(res<1e-7&&old<1e-7))
        return res;
    old=res;
  }
  return 0.;
}
