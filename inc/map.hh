#include <cmath>
/*--------------------------------------------------------------------*/
// implement variable changes:
//
//  map_Linear,   [0,1] -> [a,b]
//  map_Finite,   [0,1] -> [a,b] for endpoint singularities
//  map_SemiInf,  [0,1] -> [a,+inf] (or [-inf,a] if a<0)
//
struct map {
  double a, b;//, tmax;
  virtual double operator ()(double t) = 0; // t is the new var
};

template <typename F>
struct Finite : map {
  //double a, b, tmax=4.;
  double tmax = 3.4;
  F funk;
  double operator ()(double t) {
    double dxdt, del, Q = exp(-2.*sinh(t*tmax));
    dxdt = tmax*2.*(b-a)*Q*cosh(t*tmax)/( (1.+Q)*(1.+Q) );
    del = (b-a)*Q/(1.+Q);
    double res = ( // trick to avoid cancellation errors
        funk( a+del , del )+
        funk( b-del , del ) 
        )*dxdt;

    //if ( isinf(res)||isnan(res) ) { return 0.;}
    //else 
      return res;
  }
  Finite(F _func, double _a, double _b) :
    funk(_func) { a=_a; b=_b; }
};

template <typename F>
struct SemiInf : map {
  //double a, tmax=4.;
  double tmax = 18., tmin = -3.5;
  //double tmax = 3.8;
  F funk;
  double operator ()(double t) {
    //double dxdt, del = a*exp(2.*sinh(((1.-2.*t)*tmax)));
    //dxdt = tmax*4.*del*cosh(((1.-2.*t)*tmax));
    //double tt = (1.-2.*t)*tmax;
    //double dxdt, del = a*exp( tt - exp( -tt ) );
    //dxdt = -tmax*2.*del*(1.+exp(-tt));
    //return funk(a+del,del)*fabs(dxdt);
    
    /*double dxdt, 
           del1 = a*exp(2.*sinh((+t*tmax))),
           del2 = a*exp(2.*sinh((-t*tmax)));
    dxdt = tmax*2.*cosh((t*tmax));
    double res = 
       ( funk(a+del1,del1)*fabs(del1)
       + funk(a+del2,del2)*fabs(del2)
      )*fabs(dxdt);//*/
    
    double del1 = a*exp( +t*tmax - exp(-t*tmax)),
           del2 = a*exp( +t*tmin - exp(-t*tmin));

    double res =
           funk(a+del1,del1)*fabs(del1*(1.+exp(-t*tmax)))*fabs(tmax)
         + funk(a+del2,del2)*fabs(del2*(1.+exp(-t*tmin)))*fabs(tmin)  ;//*/

    //if ( isinf(res)||isnan(res) ) { return 0.;}
    //else 
      return res;
  }
  SemiInf(F _func, double _a) :
    funk(_func) { a=_a; }
};

// helpers:
template<typename F>
Finite<F> remap(F f,double a, double b) { return Finite<F>(f,a,b); }
template<typename F>
SemiInf<F> remap(F f,double a) { return SemiInf<F>(f,a); }
/* 
 * Example of usage with lambda function, changing variables
 * from p in the range (-1,2) into x over (0,1). The function
 * is f(p)=p^2 which becomes f(p(x))*dp/dx.
 *
 *    auto f = [&](double p, double pd) { return p*p; };
 *    remap(f, -1., 2. )(x); 
 *
 *  The second argument "pd" is not used here. It can cater for
 *  integrable singularities at the interval endpoints.
 */
