#include "minim.h"

double calc_func(double funs(int i, double x),int funsnr, vector* params, vector* errors, int sgn, double x){//sgn should be +1/-1/0 depending on the values
  double buffer=0;
  for (int i=0;i<funsnr;i++){
    buffer+=funs(i,x)*(params->data[i]+sgn*errors->data[i]);
  }
  return buffer;
}


int main(int argc, char** argv){
  int sweep;/*
  vector* calls=vector_alloc(3);
*/  printf("start with the first function -Rosenbrock's valley function- at the (0,0) point, since this is what my basic vector allocation creates\n" );
  vector* x=vector_alloc(2);
  printf("  here is the starting vector: %g \t %g \n \n", x->data[0], x->data[1]);
  sweep=newton(	 func2, x,0.000001);
  printf("solution reached in %d steps\n", sweep );
  printf("Solution is (should be (a, a^2, which in our case a=1 -> (1,1))): \n  %g %g  \n \n", x->data[0], x->data[1]);

  //printf("It took %d steps and %g functioncalls \n", sweep, calls->data[0]);
  printf("Starting the second function -Himmelblau's function-, starting point is the output of the previous search, and now solution should be: \n either x=-0.270845 y=− 0.923039 (case: local maximum) \n or one of the following:  \n 3.0 2.0 \n -2.805118 3.131312  \n -3.779310  -3.283186 \n 3.584428  -1.848126 \n" );
  sweep=newton(func3, x, 0.000001);
  printf("solution reached in %d steps\n", sweep );
//  printf("It took %d steps and %g functioncalls \n", sweep, calls->data[2]);
  printf("solution is: %g %g\n Now we found the maximum - since the gradient is zero here too \n ",  x->data[0], x->data[1] );

  printf("let's start from starting points which are near the minima - since newton steps are supposed to solve these type of problems\n" );
  printf("starting point: (2,2)\n" );
  x->data[1]=2;
  x->data[0]=2;
  sweep=newton(func3, x, 0.000001);
  printf("solution reached in %d steps\n", sweep );
//  printf("It took %d steps and %g functioncalls \n", sweep, calls->data[2]);
  printf("solution is: %g %g\n",  x->data[0], x->data[1] );
  printf("For another minimum, starting point: (-2,2)\n" );
  x->data[1]=2;
  x->data[0]=-2;
  sweep=newton(func3, x, 0.000001);
  printf("solution reached in %d steps\n", sweep );
//  printf("It took %d steps and %g functioncalls \n", sweep, calls->data[2]);
  printf("solution is: %g %g\n",  x->data[0], x->data[1] );



  printf(" \n \n I can conclude that we've succesfully found multiple extrema. \n \n Part A is finished, here comes part B\n ----------------------------------------------------------- \n " );
/////////////////////////////////////////////////////////////////////////////

/*  printf("start with the first equation system \n we will start from (0,1) point\n" );
  x->data[0]=0;
  x->data[1]=1;
  for (int i=0;i<3;i++){
    calls->data[i]=0;//reset the counter
  }
  printf("  here is the starting vector: %g \t %g \n \n", x->data[0], x->data[1]);
  sweep=newtonB(	 func1_noJ, x, 0.000001, calls, 0, 0.001);
  printf("It took %d steps and %g functioncalls \n", sweep, calls->data[0]);
  printf("solution vector (should be - with epsilon=0.00001- around 9.1061,0.00011 - order shouldn't matter) :\n" );
  printf(" %g \t %g \n \n", x->data[0], x->data[1]);
  printf("Starting the second function: We use the points above as a starting point, and should find (1,1) as solution\n" );
  sweep=newtonB(func2_noJ, x, 0.000001, calls, 1, 0.001);
  printf("It took %d steps and %g functioncalls \n", sweep, calls->data[1]);
  printf("Solution is: %g %g  \n \n", x->data[0], x->data[1]);
  printf("Starting the third function, starting point is once again the output of the previous search, and now solution should be: \n either x=-0.270845 y=− 0.923039 (case: local maximum) \n or one of the following:  \n 3.0 2.0 \n -2.805118 3.131312  \n -3.779310  -3.283186 \n 3.584428  -1.848126 \n" );

  sweep=newtonB(func3_noJ, x, 0.000001, calls, 2, 0.001);
  printf("It took %d steps and %g functioncalls \n", sweep, calls->data[2]);
  printf("solution is: %g %g\n",  x->data[0], x->data[1] );
  printf("\n We see something interesting in the functioncalls, on the 2nd case, method B seems more efficient (probably the reason is the new stop condition) for other cases, we have more functioncalls in a similar amount of steps\n" );

  vector_free(x);
*/
}
