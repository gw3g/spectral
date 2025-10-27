#pragma once
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <memory>
/*--------------------------------------------------------------------*/
// wrapper for GSL ∫ methods

inline gsl_integration_workspace* get_thread_workspace(size_t n=2000) {
    static thread_local std::unique_ptr<gsl_integration_workspace> ws{
        gsl_integration_workspace_alloc(n)
    };
    return ws.get();
}

struct quad {
  gsl_integration_workspace * wsp;

  quad(const size_t n=1000):           // recursion size
    wsp(gsl_integration_workspace_alloc(n)) {}

  ~quad() {                // free memory at destruction
    gsl_integration_workspace_free(wsp); }

  operator gsl_integration_workspace*() { return wsp; }
};

template <typename F>
class gsl_function_pp: gsl_function {
  const F func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->func(x);
  }
  public:
  gsl_function_pp(const F& f) : func(f) {
    function = &gsl_function_pp::invoke; //inherited from gsl_function
    params   = this;
  }
  operator gsl_function*(){return this;}
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
  return gsl_function_pp<F>(func);
};

