# NLO spectral functions

Here I provide functions to compute the imaginary part of 
two-loop diagrams, like the one shown below.
These are thermal self energies, functions of an external four momentum _K_
and dependent on the statistical configuration.

![Labelling of generic two-loop diagram](inc/twoloop.png?raw=true "2-loop")

Powers of the internal propagators for _P_, _Q_, _R_, _L_ and _V_ respectively are 
used to classify the masters by a string `abcde'.
Here the powers are all either 0, 1 or 2. 
(Propagators in the figure will be absent if the power is 0.)
The s<sub>i</sub> are statistical factors: +1 for bosons and -1 for fermions.
By "cutting" the diagram, some momenta are put on-shell which allows
the master to be computed as a 2-dimensional integral over _p_ and _q_
(corresponding to the 3-momenta of the loop variables above).

See below for a list of the included functions.
(They are defined in `inc/core.hh`)

## Usage

By default, the program handles a master integral that is read from the
`config` file in the working directory. 
An example, to demonstrate the correct formatting, consider:
```
11020
++-
1
0
```
This file sets up a master integral with no propagators for _R_ and _V_.
Those with _P_ and _Q_ are single propagators and _L_ is repeated.
Here s<sub>0</sub>=s<sub>1</sub>=+1 and s<sub>2</sub>=-1.
Here `m=1` and `n=0` are  powers of energies 
_p_<sub>0</sub> and _q_<sub>0</sub> respectively.

The project can be built using `make' in this directory,
which puts the executable in `bin/'.
It can then be run, for example
```
./bin/rho -k 1.0
```
will set k/T=1.0 and then perform a sweep of energies from
k<sub>0</sub>/T=10<sup>-2</sup> to 10<sup>+2</sup> using
500 points (this can be changed by modifying `main.cpp' and 
recompiling).

Use the flag `-h' for more details.
Not that any changes to `main.cpp` will require a clean to be picked up
by the makefile.

### How it works

You can use the master integrals by first declaring
a pointer to a `Master` object, and then constructing
the type with energy power and exponents:
```
Master *rho;
rho = _10120(m,n,s);
```
The `s` pointer is to an `int [3]` array of entries +1 (boson)
or -1 (fermion) to capture the statistical 
configuration: (s<sub>0</sub>,s<sub>1</sub>,s<sub>2</sub>).

To get the value of the function, one can simply call:
```
(*rho)(k0,k)
```
Some examples are provided.


**NB** Instances of the masters return a quantity in units
of the temperature.


## Data files

Output from the main function is stored under `out/data`. 
(If the file is empty, it will be created by the makefile.)
Naming convention is
```
diag.[type]{k=[k/T]}.([s0][s1][s2]).[m][n].dat
```
where `[...]` takes the value specified within the square brackets.
Each file tabulates the energy master integrals, as a function of energy.
There are four columns: 1) energy k<sub>0</sub>/T, 2) rho, scaled so that ...
3) is the leading OPE and 4) the next-to-leading OPE.

These can be plotted with the helpt of the [gnuplot](https://www.gnuplot.info)
scripts in the `out` directory.


## List of included integrals

| Type        | I            | II                         | III   | IV    | V           | VI    |
|:------------|:-------------|:---------------------------|:------|:------|:------------|:------|
| Propagatars | 01020, 00120 | 11010, 10110, 11020, 10120 | 11011 | 11100 | 11110, Star | 11111 |

