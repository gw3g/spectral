# NLO spectral functions

Here I provide functions to compute the imaginary part of 
two-loop diagrams, like the one shown below.
These are thermal self energies, functions of an external four momentum K
and dependent on the statistical configuration.

![Labelling of generic two-loop diagram](inc/twoloop.png?raw=true "2-loop")

Powers of the internal propagators for P, Q, R, L and V respectively are 
used to classify the masters by a string `abcde'.
Here the powers are all either 0, 1 or 2. 
(Some propagators in the figure will be absent if the power is 0.)
The s<sub>i</sub> are statistical factors +1 for bosons
and -1 for fermions.
By cutting the diagram, some momenta are put on shall which allows
the master to be computed as a 2D integral over p and q (3-momenta magnitudes).

See below for a list of which functions have been included.
(They are defined in `inc/core.hh`)


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
configuration: (s<sub>0</sub>,s<sub>1</sub>,s<sub>2</sub>).

To get the value of the function, one can simply call:
```
(*rho)(k0,k)
```
Any changes to `main.cpp` will require a clean to be picked up
by the makefile.

**NB** Instances of the masters return a quantity in units
of the temperature. 

### List of included integrals

| Type     |  Propagators |
|:---------|:------------:|
| I        | 01020        |
|          | 00120        |
| II       | 11010        |
|          | 10110        |
|          | 11020        |
|          | 10120        |
| III      | 11011        |
| IV       | 11100        |
| V        | 11100        |
|          | (+"Star")    |
| VI       | 11111        |


### Dileptons

See files under `spike`, compile with `make dileptons`.


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

