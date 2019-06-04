#include "root.h"

double calc_func(double funs(int i, double x),int funsnr, vector* params, vector* errors, int sgn, double x){//sgn should be +1/-1/0 depending on the values
  double buffer=0;
  for (int i=0;i<funsnr;i++){
    buffer+=funs(i,x)*(params->data[i]+sgn*errors->data[i]);
  }
  return buffer;
}


int main(int argc, char** argv){
  printf("start with the first equation system \n we will start from (0,1) point\n" );
  vector* x=vector_alloc(2);
  x->data[1]=1;
  printf("  here is the starting vector: %g \t %g \n \n", x->data[0], x->data[1]);
  newton_with_jacobian(	 func1, x, 0.000001);
  printf("solution vector (should be - with epsilon=0.00001- around 9.1061,0.00011 - order shouldn't matter) :\n" );
  printf(" %g \t %g \n \n", x->data[0], x->data[1]);
  printf("Starting the second function: We use the points above as a starting point, and should find (1,1) as solution\n" );
  newton_with_jacobian(func2, x, 0.000001);
  printf("Solution is: %g %g  \n \n", x->data[0], x->data[1]);
  printf("Starting the third function, starting point is once again the output of the previous search, and now solution should be: \n either x=-0.270845 y=− 0.923039 (case: local maximum) \n or one of the following:  \n 3.0 2.0 \n -2.805118 3.131312  \n -3.779310  -3.283186 \n 3.584428  -1.848126 \n" );
/*  x->data[1]=-3.584428;
  x->data[0]=-1.848126;*/
  newton_with_jacobian(func3, x, 0.000001);
  printf("solution is: %g %g\n",  x->data[0], x->data[1] );
  vector_free(x);

}
