/*
 Solves Laplace equation in 2d
        -Laplace(u) = f in (0,1)x(0,1)
                 u  = 0 on boundary
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>

typedef Eigen::MatrixXd Array2D;

using namespace std;

// RHS function f(x,y) in Poisson eqn
double rhs(const double x, const double y)
{
   return 1.0;
}

// Computes L2 norm of residual
double residual (const double h, const Array2D& u, const Array2D& b)
{
   const unsigned int nx = u.rows();
   const unsigned int ny = u.cols();
   const double ihh = 1.0/(h*h);
   
   double res = 0.0;
   for(unsigned int i=1; i<nx-1; ++i)
      for(unsigned int j=1; j<ny-1; ++j)
      {
         double res1 = - b(i,j)
                       + ihh * ( 4.0 * u(i,j) - u(i-1,j) - u(i+1,j) 
                                              - u(i,j-1) - u(i,j+1) );
         res += res1 * res1;
      }
   
   return sqrt(res/((nx-2)*(ny-2)));
}

// Main function
int main()
{
   const double TOL = 1.0e-6;
   const unsigned int n = 50;
   const double h = 1.0/(n-1);
   
   Array2D u(n,n), uold(n,n), b(n,n);
   
   // compute rhs for interior grid points
   for(unsigned int i=1; i<n-1; ++i)
      for(unsigned int j=1; j<n-1; ++j)
      {
         double x = i * h;
         double y = j * h;
         b(i,j) = rhs(x, y);
      }
   
   // Initial guess is zero
   u.setZero();

   const double res0 = residual (h, u, b);
   double res  = res0;
   unsigned int iter = 0;
   while (res > TOL*res0)
   {
      uold = u;
      
      for(unsigned int i=1; i<n-1; ++i)
         for(unsigned int j=1; j<n-1; ++j)
            u(i,j) =   0.25 * (uold(i-1,j) + uold(i+1,j) + 
                               uold(i,j-1) + uold(i,j+1))
                     + 0.25 * h * h * b(i,j);
      
      res = residual (h, u, b);
      ++iter;
      cout << iter << "  " << res << "  " << res/res0 << endl;
   }

   // save for gnuplot
   ofstream sol("u.dat");
   for(unsigned int i=0; i<n; ++i)
   {
      for(unsigned int j=0; j<n; ++j)
      {
         double x = i * h;
         double y = j * h;
         sol << x << "  " << y << "  " << u(i,j) << endl;
      }
      sol << endl;
   }
   cout << "Saved solution into u.dat\n";
}
