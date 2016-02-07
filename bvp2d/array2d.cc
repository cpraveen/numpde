#include "array2d.h"

using namespace std;

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
