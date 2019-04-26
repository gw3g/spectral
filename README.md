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

## Data files

Output from main functions is stored under `out/data`.
Naming convention is
```
diag.[type]{k=[k/T]}.([s0][s1][s2]).[m][n].dat
```
where `[...]` takes the value specified within the square brackets.

## Todo

* deal w/ divergence by expanding log.
* ~~type 5 + star~~ tune integrator
* ~~include OPE in masters~~ check the T^4
* organise include files better
* cuda implementation for high accuracy?
* check m,n>0 results

Remarks: should I be using Finite(F &func, a, b) in map.hh?
