#include <iostream>
#include "array2d.h"

using namespace std;

// Empty constructor
Array2D::Array2D ()
:
nrow (0),
ncol (0),
n  (0)
{
}

// Constructor based on size
Array2D::Array2D (const unsigned int nrow, const unsigned int ncol)
:
nrow (nrow),
ncol (ncol),
n  (nrow*ncol),
u  (nrow*ncol)
{
#if defined(DEBUG)
   if(nrow < 0 || ncol < 0)
   {
      cout << "Dimensions cannot be negative" << endl;
      exit(0);
   }
#endif
}

// Change size of array
void Array2D::resize(const unsigned int nrow1, const unsigned int ncol1)
{
#if defined(DEBUG)
   if(nrow1 < 0 || ncol1 < 0)
   {
      cout << "Dimensions cannot be negative" << endl;
      exit(0);
   }
#endif
   nrow = nrow1;
   ncol = ncol1;
   n  = nrow1 * ncol1;
   u.resize (n);
}

// return number of rows, size of first index
unsigned int Array2D::rows() const
{
   return nrow;
}

// return number of columns, size of second index
unsigned int Array2D::cols() const
{
   return ncol;
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
   if(nrow != a.rows() || ncol != a.cols())
   {
      cout << "Array sizes do not match" << endl;
      exit(0);
   }
#endif
   u = a.u;
   return *this;
}
