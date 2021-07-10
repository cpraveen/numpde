#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>

class Array2D
{
public:
   Array2D();
   Array2D(const unsigned int nx, const unsigned int ny);
   void resize (const unsigned int nx, const unsigned int ny);
   unsigned int sizex() const;
   unsigned int sizey() const;

   inline const double&  operator()(const unsigned int i, const unsigned int j) const
   {
      #if defined(DEBUG)
      if (i < 0 || i > nx - 1 || j < 0 || j > ny - 1)
      {
         cout << "Indices out of range" << endl;
         cout << "i, j = " << i << ", " << j << endl;
         exit(0);
      }
      #endif
      return u[i * ny + j];
   }

   inline double& operator()(const unsigned int i, const unsigned int j)
   {
      #if defined(DEBUG)
      if (i < 0 || i > nx - 1 || j < 0 || j > ny - 1)
      {
         cout << "Indices out of range" << endl;
         cout << "i, j = " << i << ", " << j << endl;
         exit(0);
      }
      #endif
      return u[i * ny + j];
   }

   Array2D& operator= (const double scalar);
   Array2D& operator= (const Array2D& u);

private:
   unsigned int nx, ny, n;
   std::vector<double> u;
};

#endif
