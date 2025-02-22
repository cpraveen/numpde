/*
 * Solve 
 *      -u'' = (2*pi)^2 sin(2*pi*x) in (0,1)
 *      u(0) = u0, u(1) = u1
 * Exact solution
 *       u(x) = sin(2*pi*x) + (1-x)*u0 + x*u1
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "Vector.h"
#include "sparse_matrix.h"
#include "cg_solver.h"
#include "jacobi_solver.h"
#include "sor_solver.h"
#include "ssor_solver.h"
#include "tridiag_solver.h"

using namespace std;

//------------------------------------------------------------------------------
// Problem definition and exact solution
//------------------------------------------------------------------------------
double u0 = 0, u1 = 1;
double xmin = 0, xmax = 1;
double exact_solution(const double x)
{
   return sin(2*M_PI*x) + (1-x)*u0 + x*u1;
}

//------------------------------------------------------------------------------
// Main program
//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
   if(argc <= 2)
   {
      cout << "Specify: n, solver (jacobi, sor, ssor, cg), max_iter\n";
      cout << "Example: " << argv[0] << " 100 jacobi 1000\n";
      exit(1);
   }

   const unsigned int n = atoi(argv[1]);
   const string method = string(argv[2]);
   unsigned int max_iter = 1000;
   if(argc == 4) max_iter = atoi(argv[3]);

   // Create 1-d mesh
   const double h = (xmax - xmin) / (n - 1);
   Vector<double> x(n);
   for(unsigned int i=0; i<n; ++i)
      x(i) = xmin + i*h;

   const double a =  2.0/(h*h);
   const double b = -1.0/(h*h);

   // Construct matrix
   SparseMatrix<double> A(n);
   A.set(0, 0, a); 
   for(unsigned int i=1; i<n-1; ++i)
   {
      A.set(i, i  , a); // diagonal element should be first
      A.set(i, i-1, b);
      A.set(i, i+1, b);
   }
   A.set(n-1, n-1, a);
   A.close();

   // Construct right hand side vector
   Vector<double> f(n);
   for(unsigned int i=0; i<n; ++i)
      f(i) = pow(2*M_PI,2) * sin(2*M_PI*x(i));

   // Apply boundary conditions
   f(0)   = A(0,0) * u0;
   f(n-1) = A(n-1,n-1) * u1;

   // make symmetric
   f(1)      -= A(1,0) * u0;
   f(n-2)    -= A(n-2,n-1) * u1;
   A(1,0)     = 0;
   A(n-2,n-1) = 0;

   // Print for debugging
   //cout << "A = \n" << A << endl;
   //cout << "f = \n" << f << endl;

   // Solution vector
   Vector<double> u(n);
   u      = 0;
   u(0)   = u0;
   u(n-1) = u1;

   // Create solver object
   const double tol = 1.0e-6;
   const double omega = 2.0/(1.0 + sin(M_PI*h)); // SOR relaxation

   unsigned int iter = 0;
   if(method == "jacobi")
   {
      JacobiSolver<double> solver (max_iter, tol);
      iter = solver.solve (A, u, f);
   }
   else if(method == "sor")
   {
      SORSolver<double> solver (max_iter, tol, omega);
      iter = solver.solve (A, u, f);
   }
   else if(method == "ssor")
   {
      SSORSolver<double> solver (max_iter, tol, omega);
      iter = solver.solve (A, u, f);
   }
   else if(method == "cg")
   {
      CGSolver<double> solver (max_iter, tol);
      iter = solver.solve (A, u, f);
   }
   else if(method == "tri")
   {
      TriDiagSolver<double> solver;
      iter = solver.solve (A, u, f);
   }
   else
   {
      cout << "Unknown solver: " << method << "\n";
      exit(1);
   }
   
   cout << "Convergence tolerance = " << tol << endl;
   cout << "Number of iterations = " << iter << endl;

   // Save solution to file
   string fname = "sol.dat";
   ofstream fsol(fname);
   for(unsigned int i=0; i<n; ++i)
      fsol << x(i) << "  " 
           << u(i) << "  " 
           << exact_solution(x(i)) << endl;
   fsol.close ();
   cout << "Saved solution into file " << fname << endl;
}
