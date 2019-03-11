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
        funk( a+del , del )+
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
    double dxdt, del = a*exp(2.*sinh(((1.-2.*t)*tmax)));
    dxdt = tmax*4.*del*cosh(((1.-2.*t)*tmax));
    return funk(a+del,del)*dxdt;
  }
  SemiInf(F &_func, double _a) :
    funk(_func) { a=_a; }
};

