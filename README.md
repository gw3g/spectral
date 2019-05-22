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
configuration: (s<sub>0</sub>,s<sub>1</sub>,<sub>s</sub>).

To get the value of the function, one can simply call:
```
(*rho)(k0,k)
```

## Data files

Output from the main function is stored under `out/data`. 
(If the file is empty, it will be created by the makefile.)
Naming convention is
```
diag.[type]{k=[k/T]}.([s0][s1][s2]).[m][n].dat
```
where `[...]` takes the value specified within the square brackets.

## Todo

* type 6, m=2
* ~~type 5 + star~~ tune integrator
* ~~include OPE in masters~~ check the T^4
* ~~organise include files better~~
* cuda implementation for high accuracy?

