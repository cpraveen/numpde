#include <iostream>
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
#if defined(DEBUG)
   if(nx < 0 || ny < 0)
   {
      cout << "Dimensions cannot be negative" << endl;
      exit(0);
   }
#endif
}

// Change size of array
void Array2D::resize(const unsigned int nx1, const unsigned int ny1)
{
#if defined(DEBUG)
   if(nx1 < 0 || ny1 < 0)
   {
      cout << "Dimensions cannot be negative" << endl;
      exit(0);
   }
#endif
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
#if defined(DEBUG)
   if(i < 0 || i > nx-1 || j < 0 || j > ny-1)
   {
      cout << "Indices out of range" << endl;
      cout << "i, j = " << i << ", " << j << endl;
      exit(0);
   }
#endif
   return u[i + j*nx];
}

// Return reference to (i,j), this can modify the value
double& Array2D::operator() (const unsigned int i, const unsigned int j)
{
#if defined(DEBUG)
   if(i < 0 || i > nx-1 || j < 0 || j > ny-1)
   {
      cout << "Indices out of range" << endl;
      cout << "i, j = " << i << ", " << j << endl;
      exit(0);
   }
#endif
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
#if defined(DEBUG)
   if(nx != a.sizex() || ny != a.sizey())
   {
      cout << "Array sizes do not match" << endl;
      exit(0);
   }
#endif
   u = a.u;
   return *this;
}
