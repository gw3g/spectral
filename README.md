# NLO spectral functions

## Usage

You can use the master integrals by first declaring
a pointer to a `Master` object, and then constructing
the type with energy power and exponents:
```
Master *rho;
rho = _10120(m,n,s);
```
Here `m` and `n` are the powers of p_0 and q_0 (resp.),
and `s` is an `int [3]` array of entries +1 (boson)
or -1 (fermion) to capture the statistical 
configuration: (s_0,s_1,s_2).

To get the value of the function, one can simply call:
```
(*rho)(k0,k)
```

## Todo

* type 4,5,6
* include OPE in masters
* 2D integrator {Trapezoidal,GSL Gauss-Konrad}
* organise include files better
* cuda implementation for high accuracy?

Remarks: should I be using Finite(F &func, a, b) in map.hh?
