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
// runner, eps is relative
#include <iostream>
using namespace std;

template <typename F>
double go(integrate<F> &I, double eps=1.e-4) {
  const int iMax=10; // max number of iterations
  double res, old;
  for (int i=0;i<iMax;i++) {
    res = I.next(); //cout << "Step: " << i << "  \
          Result = "<< res << endl;
    if (i>5)
      if (abs(res-old)<eps*abs(old)||(res<1e-7&&old<1e-7))
        return res;
    old=res;
  }
  return 0.;
}
/*
// for 2D integration
template <typename F>
struct inner {
  double _x;
  F *func;
  double operator ()(double y) { return (*func)(_x,y); }
};

template <typename F>
struct outer {
  inner<F> f2;
  double operator ()(double x) {
    f2._x = x;
    integrate<inner<F>> I(f2);
    return go(I);
  }
};

template <typename F>
double integrate_2d(F *_func) {
  outer<F> f1;
  f1.f2.func = _func;
  integrate<outer<F>> I(f1);
  return go(I);
}*/
