#ifndef __TRIDIAG_SOLVER_H__
#define __TRIDIAG_SOLVER_H__

#include "sparse_matrix.h"
#include "Vector.h"

template <class T>
class TriDiagSolver
{
   public:
      TriDiagSolver ();
      ~TriDiagSolver () {};
      unsigned int solve (const SparseMatrix<T>& A,
                                      Vector<T>& x, 
                          const       Vector<T>& f);

   private:
      void lu_factor();
      void lu_solve(Vector<T>& x, const Vector<T>& f) const;

      bool factored;
      Vector<double> a, b, c;
};

#endif
