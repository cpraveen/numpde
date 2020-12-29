# Linear advection in 2d

Authors: Ashish Bhole, Praveen C

There are two methods implemented: 

* FVM with upwind flux, many limiters and RK time stepping
* Multi-D upwind scheme of LeVeque
* Lax-Wendroff scheme

See the input file `run/inp.dat` for more parameters.

The method of LeVeque is from

LeVeque: High resolution conservative algorithms for advection in incompressible flows, SIAM J. NUMER. ANAL, Vol. 33, No. 2, pp. 627-665, April 1996. https://doi.org/10.1137/0733033

## How to run the code

Compile the code in `src`

```
cd src
make
```

This creates the executable `advect`

Edit input file `inp.dat` in `run` directory

Run the code inside `run` directory

```
../src/advect
```

Visualize solutions in Visit

## Convergence rate

Edit the python code `conv.py` to set scheme parameters

Run the python code

```
python conv.py
```

It plots error versus N and also saves the plot into `error.pdf`
