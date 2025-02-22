#include <cassert>
#include <cmath>

#include "tridiag_solver.h"
#include "math_functions.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
TriDiagSolver<T>::TriDiagSolver ()
{
   std::cout << "Using TriDiag solver\n";
   factored = false;
}

//-----------------------------------------------------------------------------
// Perform LU factorization
//-----------------------------------------------------------------------------
template <class T>
void TriDiagSolver<T>::lu_factor ()
{
   assert(factored == false);

   // Do the factorization

   factored = true;
}

//-----------------------------------------------------------------------------
// Solve the matrix
//-----------------------------------------------------------------------------
template <class T>
void TriDiagSolver<T>::lu_solve (Vector<T>& x, const Vector<T>& f) const
{
   assert(factored == true);
}
      
//-----------------------------------------------------------------------------
// Solves A*x = f for x
// We assume that x has already been initialized.
//-----------------------------------------------------------------------------
template <class T>
unsigned int TriDiagSolver<T>::solve (const SparseMatrix<T>& A,
                                                  Vector<T>& x,
                                      const       Vector<T>& f)
{
   const unsigned int n = x.size();
   assert (n == A.size());
   assert (n == f.size());

   // Allocate memory for a,b,c

   // Copy diagonals of A into a,b,c

   // Perform LU decomposition
   lu_factor();

   // Perform backward and forward solve
   lu_solve(x, f);

   return n;
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class TriDiagSolver<float>;
template class TriDiagSolver<double>;
