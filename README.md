# NLO spectral functions

Based on: [1910.07552](https://arxiv.org/abs/1910.07552) (PRD,2019).

Here I provide functions to compute the imaginary part of 
two-loop diagrams, like the one shown below.
These are thermal self energies, functions of an external four momentum 
_K_=(k<sub>0</sub>,**k**)
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

## Physical applications

Two specific examples are supplied.

### QGP dilepton rate

See: [1910.09567](https://arxiv.org/abs/1910.09567)

The combination of master integrals which is required to compute the dilpton emission rate
from a hot QCD plasma is demonstrated in `examples/dileptons.cpp` (see functions in file for more information).
This source file can be compiled with:
```
make dileptons
```
This file was used to compute the rates implemented in hydrodynamic simulations in
 [2311.06675](https://arxiv.org/abs/2311.06675) and [2311.06951](https://arxiv.org/abs/2311.06951).
Note that these spectral functions are valid for finite quark chemical potential.


### Neutrino interaction rate

See: [2312.07015](https://arxiv.org/abs/2312.07015)

The QED corrections to the neutrino interaction rate at MeV temperatures
requires an integration over the electron-positron spectral function.
We supply code for evaluating this rate in `examples/neutrinos.cpp` (see functions in file, accuracy
may be adjusted to suit).
This source file can be compiled with:
```
make neutrinos
```


## Usage

The program handles a master integral that is read from the
`config` file in the working directory. 
An example, which demonstrates the correct formatting, reads:
```
11020
++-
1
0
```
When processed, this input sets up a master integral with 
no propagators for _R_ and _V_.
Those with _P_ and _Q_ are single propagators and _L_ is repeated.
Here s<sub>0</sub>=s<sub>1</sub>=+1 and s<sub>2</sub>=-1 summarise
the statistical
configuration: (s<sub>0</sub>,s<sub>1</sub>,s<sub>2</sub>).
Here `m=1` and `n=0` are  powers of energies 
_p_<sub>0</sub> and _q_<sub>0</sub> respectively.

The main project can be built using `make` in the same directory as 
this readme, which puts the executable in `bin/`.
It can then be run, for example using
```
./bin/rho -k 1.0
```
which will set k/T=1.0 and then perform a sweep of energies from
k<sub>0</sub>/T=10<sup>-2</sup> to 10<sup>+2</sup> using
500 points on a logarithmic scale 
(this can be changed by modifying `main.cpp` and recompiling).

Use the flag `-h` to see more options.
Note that any changes to `main.cpp` will require a clean to be picked up
by the makefile.

### Data files

Output from the main function is stored under `out/data`. 
(If the file is empty, it will be created by the makefile.)
Naming convention is
```
diag.[type]{k=[k/T]}.([s0][s1][s2]).[m][n].dat
```
where `[...]` takes the value specified within the square brackets.
Each file tabulates the energy master integrals, as a function of energy.
There are four columns: (1) energy k<sub>0</sub>/T, (2) rho scaled with T=1,
(3) the leading OPE and (4) the next-to-leading OPE.

These can be plotted with the help of the [gnuplot](https://www.gnuplot.info)
scripts in the `out` directory.


## Abusage

You can manipulate the master integrals by first declaring
a pointer to a `Master` object, and then constructing
the type with energy power and exponents:
```
Master *rho;
rho = _11020(m,n,s);
```
The `s` pointer is to an `int [3]` array of entries +1 (boson)
or -1 (fermion).

To get the value of the function, one can simply call:
```
(*rho)(k0,k)
```
Some examples are provided.


**NB**, Instances of the master functions return a quantity in units
of the temperature (i.e. T=1).


## Requirements

* [GNU Science Library](https://www.gnu.org/software/gsl)

* C++ compiler (v11 or higher)



## List of included integrals

| Type        | I            | II                         | III   | IV    | V           | VI    |
|:------------|:-------------|:---------------------------|:------|:------|:------------|:------|
| Propagators | 01020, 00120 | 11010, 10110, 11020, 10120 | 11011 | 11100 | 11110, Star | 11111 |

