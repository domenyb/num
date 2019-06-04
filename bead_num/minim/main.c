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

  printf("start with the first function -Rosenbrock's valley function- at the (0,0) point, since this is what my basic vector allocation creates\n" );
  x->data[0]=0;
  x->data[1]=0;
  printf("  here is the starting vector: %g \t %g \n \n", x->data[0], x->data[1]);
  sweep=qsi_newton(	 func2_noH, x,0.000001);
  printf("solution reached in %d steps\n", sweep );
  printf("Solution is (should be (a, a^2, which in our case a=1 -> (1,1))): \n  %g %g  \n \n", x->data[0], x->data[1]);
  vector_free(x);
  printf("I can't seem to find the problem - most likely a missed + or - in an equation, and I don't have enough time left to find it, so altough the excersise is almost done, I sadly have to settle with minimal result\n" );

  //printf("It took %d steps and %g functioncalls \n", sweep, calls->data[0]);
/*  double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
  double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
  double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
  int N = sizeof(t)/sizeof(t[0]);

  vector* x2=vector_alloc(3);
  x2->data[0]=3;
  x2->data[1]=1;
  x2->data[2]=0;
  printf("  here is the starting vector: %g \t %g \t %g \n \n", x2->data[0], x2->data[1], x2->data[2]);
//  sweep=qsi_newton(custom_func, x2,0.000001);
  printf("solution reached in %d steps\n", sweep );
  printf("Solution is: \n  %g \t %g \t %g  \n \n", x2->data[0], x2->data[1], x2->data[3]);
  FILE* f=fopen("data.dat", "w");
  for (double i=0;i<10; i+=0.01){
    double buffer=x2->data[0]*exp(-i/x2->data[1])+x2->data[2];
    fprintf(f, "%g \t %g \n", i, buffer );
  }
  printf("solution vs. experimental data is ready to fit \n");
  fclose(f);
  vector_free(x2);*/
}
