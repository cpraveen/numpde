#include "array2d.h"

using namespace std;

// Empty constructor
Array2D::Array2D ()
:
nx (0),
ny (0),
n  (0)
{
}

// Constructor based on size
Array2D::Array2D (const unsigned int nx, const unsigned int ny)
:
nx (nx),
ny (ny),
n  (nx*ny),
u  (nx*ny)
{
}

// Change size of array
void Array2D::resize(const unsigned int nx1, const unsigned int ny1)
{
   nx = nx1;
   ny = ny1;
   n  = nx1 * ny1;
   u.resize (n);
}

// return number of rows, size of first index
unsigned int Array2D::sizex() const
{
   return nx;
}

// return number of columns, size of second index
unsigned int Array2D::sizey() const
{
   return ny;
}

// Return value at (i,j), this is read only
double Array2D::operator() (const unsigned int i, const unsigned int j) const
{
   return u[i + j*nx];
}

// Return reference to (i,j), this can modify the value
double& Array2D::operator() (const unsigned int i, const unsigned int j)
{
   return u[i + j*nx];
}

// Set all elements to scalar value
Array2D& Array2D::operator= (const double scalar)
{
   for (unsigned int i=0; i<n; ++i)
      u[i] = scalar;
   return *this;
}

// Copy array a into this one
Array2D& Array2D::operator= (const Array2D& a)
{
   u = a.u;
   return *this;
}
