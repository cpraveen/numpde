#include <iostream>

#include "sparse_matrix.h"
#include "Vector.h"

using namespace std;

int main ()
{
   unsigned int nrow=4;
   SparseMatrix<double> A(nrow);

   A.set(0,0,10); // diagonal
   A.set(0,3,7);

   A.set(1,1,1); // diagonal
   A.set(1,2,2);

   A.set(2,2,5); // diagonal
   A.set(2,0,3);
   A.set(2,3,9);

   A.set(3,3,1); // diagonal

   A.close();

   cout << A << endl;
   cout << "A(2,3) = " << A(2,3) << endl << endl;
   
   Vector<double> x(nrow), y(nrow);

   x = 1.0;
   cout << "x = \n" << x << endl;

   A.multiply(x, y);
   cout << "y = \n" << y;
}
