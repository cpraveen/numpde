#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

class Array2D
{
public:
   Array2D();
   Array2D(const unsigned int nx, const unsigned int ny);
   void resize (const unsigned nx, const unsigned int ny);
   unsigned int sizex() const;
   unsigned int sizey() const;
   double  operator()(const unsigned int i, const unsigned int j) const;
   double& operator()(const unsigned int i, const unsigned int j);
   Array2D& operator= (const double scalar);
   Array2D& operator= (const Array2D& u);
   
private:
   unsigned int nx, ny;
   vector<double> u;
};

Array2D::Array2D (const unsigned int nx, const unsigned int ny)
:
nx (nx),
ny (ny)
{
   u.resize (nx*ny);
}

void Array2D::resize(const unsigned int nx, const unsigned int ny)
{
   u.resize (nx*ny);
}

unsigned int Array2D::sizex() const
{
   return nx;
}

unsigned int Array2D::sizey() const
{
   return ny;
}

double Array2D::operator() (const unsigned int i, const unsigned int j) const
{
   return u[i + j*nx];
}

double& Array2D::operator() (const unsigned int i, const unsigned int j)
{
   return u[i + j*nx];
}

Array2D& Array2D::operator= (const double scalar)
{
   for (unsigned int i=0; i<nx*ny; ++i)
      u[i] = scalar;
   return *this;
}

Array2D& Array2D::operator= (const Array2D& a)
{
   for (unsigned int i=0; i<nx*ny; ++i)
      u[i] = a.u[i];
   return *this;
}

double rhs(const double x, const double y)
{
   return 1.0;
}

double residual (const double h, const Array2D& u, const Array2D& b)
{
   const unsigned int nx = u.sizex();
   const unsigned int ny = u.sizey();
   const double ihh = 1.0/(h*h);
   
   double res = 0.0;
   for(unsigned int i=1; i<nx-1; ++i)
      for(unsigned int j=1; j<ny-1; ++j)
      {
         double res1 = - b(i,j) + ihh * ( 4.0 * u(i,j) - u(i-1,j) - u(i+1,j) - u(i,j-1) - u(i,j+1) );
         res += res1 * res1;
      }
   
   return sqrt(res/((nx-2)*(ny-2)));
}

int main()
{
   const unsigned int n = 50;
   const double h = 1.0/(n-1);
   
   Array2D u(n,n), uold(n,n), b(n,n);
   
   // compute rhs for interior grid points
   for(unsigned int i=1; i<n-1; ++i)
      for(unsigned int j=1; j<n-1; ++j)
      {
         double x = i*h;
         double y = j*h;
         b(i,j) = rhs(x, y);
      }
   
   u = 0.0;
   double res = residual (h, u, b);
   unsigned int iter = 0;
   while (res > 1.0e-6)
   {
      uold = u;
      
      for(unsigned int i=1; i<n-1; ++i)
         for(unsigned int j=1; j<n-1; ++j)
            u(i,j) = 0.25 * (uold(i-1,j) + uold(i+1,j) + uold(i,j-1) + uold(i,j+1))
                     + 0.25 * h * h * b(i,j);
      
      res = residual (h, u, b);
      ++iter;
      cout << iter << "  " << res << endl;
   }

   // save for gnuplot
   ofstream sol("u.dat");
   for(unsigned int i=0; i<n; ++i)
   {
      for(unsigned int j=0; j<n; ++j)
      {
         double x = i*h;
         double y = j*h;
         sol << x << "  " << y << "  " << u(i,j) << endl;
      }
      sol << endl;
   }

}