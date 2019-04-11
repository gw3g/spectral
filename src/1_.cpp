#include "core.hh"

struct rhoTrivial : Master {
  // no dependence on K -> no discontinuity
  // (this function demonstrates the inheritance structure)
 double eval() { return 0.; }
 rhoTrivial(int _m, int _n, int _s[3]) : Master(_m,_n,_s) { type=1; }
};
// function for MAIN
Master* _01020(int m, int n, int s[3]) {
  Master *R =  new rhoTrivial(m,n,s); return R;
}
Master* _00120(int m, int n, int s[3]) {
  Master *R =  new rhoTrivial(m,n,s); return R;
}

