#include "core.hh"

struct rhoTrivial : master {
  // no dependence on K -> no discontinuity
  // (this function demonstrates the inheritance structure)
 double eval() { return 0.; }
 rhoTrivial(int _m, int _n, int _s[3]) : master(_m,_n,_s) {}
};
// function for MAIN
master* _01020(int m, int n, int s[3]) {
  master *R =  new rhoTrivial(m,n,s); return R;
}
master* _00120(int m, int n, int s[3]) {
  master *R =  new rhoTrivial(m,n,s); return R;
}

