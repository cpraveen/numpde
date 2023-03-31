#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>

class Array2D
{
public:
   Array2D();
   Array2D(const unsigned int nrow, const unsigned int ncol);
   void resize (const unsigned int nrow, const unsigned int ncol);
   unsigned int rows() const;
   unsigned int cols() const;

   inline const double&  operator()(const unsigned int i, 
                                    const unsigned int j) const
   {
      #if defined(DEBUG)
      if (i < 0 || i > nrow - 1 || j < 0 || j > ncol - 1)
      {
         cout << "Indices out of range" << endl;
         cout << "i, j = " << i << ", " << j << endl;
         exit(0);
      }
      #endif
      return u[i * ncol + j];
   }

   inline double& operator()(const unsigned int i, 
                             const unsigned int j)
   {
      #if defined(DEBUG)
      if (i < 0 || i > nrow - 1 || j < 0 || j > ncol - 1)
      {
         cout << "Indices out of range" << endl;
         cout << "i, j = " << i << ", " << j << endl;
         exit(0);
      }
      #endif
      return u[i * ncol + j];
   }

   Array2D& operator= (const double scalar);
   Array2D& operator= (const Array2D& u);

private:
   unsigned int nrow, ncol, n;
   std::vector<double> u;
};

#endif
