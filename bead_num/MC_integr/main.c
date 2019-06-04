#include "integration.h"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_blas.h>
#define PI 3.141592653589793238462
double my_function(double x){return x*x;}//according to the excersise, we should use  nested functions

int main(int argc, char** argv){//most of this function is shamelessly copied from my previous excersise and is just modified for this one
  printf("A, we'd have to test this on some 'interesting integrals', well, might as well use the ones from the previous excersise\n" );
  double f1(vector* x){return sqrt(x->data[0]);};//modified to take the vector form
  int points = 1000000;
  int dim=1;
  double err=0;//we would need this to be added to be able to get it out
  vector* a=vector_alloc(dim);
  vector* b=vector_alloc(dim);
  double sol1;
  a->data[0]=0;
  b->data[0]=1;//the borders for the integration - yea, 1d, very original...
  plainmc( f1,a, b, points, &sol1, &err);
  printf("First integral from 0-1 sqrt(x) supposed to be = 2/3\n" );
  printf("The integral's value (after using 10^6 points) is: %g +/- %g \n",sol1, err );


  double f2(vector* x){return 1.0/sqrt(x->data[0]);};
  points = 1000000;
   dim=1;
   err=0;//we would need this to be added to be able to get it out
  double sol2;
  a->data[0]=0;
  b->data[0]=1;//the borders for the integration - yea, 1d, very original...
  plainmc( f2,a, b, points, &sol2, &err);
  printf("Second integral from 0-1 1/sqrt(x) is supposed to be = 2\n" );
  printf("The integral's value (after using 10^6 points) is: %g +/- %g \n",sol2, err );




  double f3(vector* x){return log(x->data[0])/sqrt(x->data[0]);};
  points = 1000000;
   dim=1;
   err=0;//we would need this to be added to be able to get it out
  double sol3;
  a->data[0]=0;
  b->data[0]=1;//the borders for the integration - yea, 1d, very original...
  plainmc( f3,a, b, points, &sol3, &err);
  printf("third integral from 0-1 ln(x)/sqrt(x) supposed to be = -4\n" );
  printf("The integral's value (after using 10^6 points) is: %g +/- %g \n",sol3, err );

  printf("Last integral, designed for this calc: ∫0π  dx/π ∫0π  dy/π ∫0π  dz/π [1-cos(x)cos(y)cos(z)]-1 = Γ(1/4)4/(4π3) ≈ 1.3932039296856768591842462603255\n" );


  double f4(vector* x){return 1.0/(PI*PI*PI*(1.0-cos(x->data[0])*cos(x->data[1])*cos(x->data[2])));};
  points = 100000000;
  dim=3;
  err=0;//we would need this to be added to be able to get it out
  double sol4;
  vector* a_big=vector_alloc(dim);
  vector* b_big=vector_alloc(dim);
  a_big->data[0]=0.0;
  a_big->data[1]=0.0;
  a_big->data[2]=0.0;
  b_big->data[0]=PI;
  b_big->data[1]=PI;
  b_big->data[2]=PI;
  plainmc(f4,a_big, b_big, points, &sol4, &err);
  printf("The integral's value (after using 10^7 points this time) is: %g +/- %g \n",sol4, err );




  printf("this is successful so far... Let's show that the Monte-Carlo method behaves as O(1/√N). \n ----------------------------------------------------------------------\n" );
  printf("Let's go back to the first integral\n" );
  FILE* f =fopen("data.dat", "w");
  for (int i=1000; i<1000000000;i*=1.78){
    dim=1;
    err=0;//we would need this to be added to be able to get it out
    a->data[0]=0;
    b->data[0]=1;//the borders for the integration - yea, 1d, very original...
    plainmc( f1,a, b, i, &sol1, &err);
    printf("First integral from 0-1 sqrt(x) supposed to be = 2/3\n" );
    printf("The integral's value (after using %d points) is: %g +/- %g \n",i, sol1, err );
    fprintf(f, "%g \t %g\n",1/sqrt(i), err );
  }


  fclose(f);
  vector_free(a);
  vector_free(b);
  vector_free(a_big);
  vector_free(b_big);

  return 0;
}
