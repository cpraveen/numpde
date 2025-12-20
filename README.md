# Numerical methods for PDE

Codes for my numerical methods for PDE course, taught at TIFR-CAM. These are based on finite difference and finite volume method.

There are some examples of boundary value problems in my [Numerical Analysis](http://github.com/cpraveen/na) codes, under numerical linear algebra.

Codes on finite element method are available [here](http://github.com/cpraveen/fem).

* fd: finite difference approximation
* bvp1d: 1D boundary value problems
* bvp2d: 2D boundary value problems
* linhyp1d: Linear hyperbolic pde in 1D
* linadv2d: linear hyperbolic pde in 2D
* claw1d: non-linear conservations laws in 1D
* heat1d: Heat equation in 1D
* euler1d: compressible Euler eqns in 1D
* euler2d: compressible Euler eqns in 2D
* multigrid: Multigrid method for BVP
* vte2d: vorticity-streamfunction form of Navier-Stokes
* matlab: older matlab examples

## Jupyter notebooks

You can copy and paste the python code into a jupyter notebook, e.g., on [colab](http://colab.research.google.com). Add following lines at the top of each notebook to get better images

```
%config InlineBackend.figure_format = 'svg'
```

## Finite Element Method

See these other repositories for FEM code

* [FEM codes in Python/deal.II/Firedrake/Fenics](https://github.com/cpraveen/fem)
* [Matlab FEM code](https://github.com/cpraveen/fem50)
* [Julia FEM code](https://github.com/cpraveen/juliafem)

---

* `Origin`: https://codeberg.org/cpraveen/numpde
* `Mirror`: https://git.sr.ht/~cpraveen/numpde
* `Mirror`: https://github.com/cpraveen/numpde