#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>

using namespace std;

class Array2D
{
public:
   Array2D();
   Array2D(const unsigned int nx, const unsigned int ny);
   void resize (const unsigned int nx, const unsigned int ny);
   unsigned int sizex() const;
   unsigned int sizey() const;
   double  operator()(const unsigned int i, const unsigned int j) const;
   double& operator()(const unsigned int i, const unsigned int j);
   Array2D& operator= (const double scalar);
   Array2D& operator= (const Array2D& u);
   
private:
   unsigned int nx, ny, n;
   vector<double> u;
};

#endif