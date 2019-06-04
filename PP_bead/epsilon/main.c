#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "extra.h"
int max_int=INT_MAX;
int small_int=INT_MIN;
float flt = FLT_EPSILON;
double dbl = DBL_EPSILON;
long double ld = LDBL_EPSILON;



int main(){
   printf("this is the largest integer %d\n", max_int);
   printf("smallest: %d\n", small_int);
   int i=1;
   while (i+1>i) {
     i++;
//     printf("%d \n", i );
   }
   printf("This is what the while reaches: %d\n",i);
    i = 0;
   while (i-1<i) {i--;}
   printf("and down: = %d\n",i);
  for (i = 0; i < i+1; i++) {}
  printf("this is what the for reaches %d\n",i);
  for (i = 0; i-1 < i; i--) {}
  printf("down: %d\n",i);

  i = 0;
  do {i++;} while(i < i+1);
  printf("dowhile: %d\n",i);
  i = 0;
  do {i--;} while(i-1 < i);
  printf("down: %d\n",i);
  float x=1; while(1+x!=1){
    x/=2;
  }
    x*=2;
  float x2=1;
  for (x2 = 1; 1+x2!=1; x2/=2){}
  x2*=2;
  float x3=1;
   do {x3/=2;} while(1+x3!=1);
   x3*=2;
  printf("float epsilon: %g\n", flt);
  printf("with while: %g\n", x);
  printf("for: %g\n", x2);
  printf("dowhile: %g\n\n", x3);

  double y1=1;double y2=1;double y3=1;
   for (y2= 1; 1+y2!=1; y2/=2){}
   y2*=2;
   while(1+y1!=1){y1/=2;}
   y1*=2;
   do {y3/=2;} while(1+y3!=1);
   y3*=2;
  printf("double epsilon: %g\n", dbl);
  printf("with while: %g\n", y1);
  printf("for: %g\n", y2);
  printf("dowhile: %g\n\n", y3);



  long double z1=1;long double z2=1;long double z3=1;
  while(1+z1!=1){z1/=2;} z1*=2;
  for (z2 = 1; 1+z2!=1; z2/=2){} ;
  z2*=2;
  do {z3/=2;} while(1+z3!=1);
  z3*=2;
  printf("long double epsilon: %Lg\n", ld);
  printf("with while: %Lg\n", z1);
  printf("for: %Lg\n", z2);
  printf("dowhile: %Lg\n\n", z3);


  int max=INT_MAX/3;
  float buff=1.0;
  float sum_up_float = 0;
  float sum_down_float = 0;
  for (i = 1.0; i < max; i++) {
    sum_up_float=sum_up_float+buff/i;
  }
  for (i = 0.0; i < max; i++) {
    sum_down_float=sum_down_float+buff/(max-i);
  }
  printf("this is the max/3: %d\n", max);
  printf("up: %g\n", sum_up_float);
  printf("down:  %g\n", sum_down_float);

  printf("The difference is  that the down-sum will add small numbers first, and after a long time it will 'disappear' \n");
  printf("BUT if we start with the big numbers, the small ones won't even matter in the first place later...\n");

  printf("With a non-infty max it should converge.\n");
  double buff2=1.0;
  double sum_up_double = 0;
  double sum_down_double = 0;
  for (i = 1.0; i < max; i++) {
    sum_up_double=sum_up_double+buff2/i;
  }
  for (i = 0.0; i < max; i++) {
    sum_down_double=sum_down_double+buff2/(max-i);
  }

  printf("Now double: \n \n Up:= %g\n", sum_up_double);
  printf("down: %g\n", sum_down_double);


  printf("test our func\n");
  double a = 1;
  double b = 2;
  double tau = 3;
  double epsilon = 12;
  int result = equal(a,b,tau,epsilon);
  printf("now we get a 1: %d\n",result);
  b = 3;
  tau = 4;
  epsilon = 1;
  result = equal(a,b,tau,epsilon);
  printf("and again: %d\n",result);

  b =23;
  tau = 1;
  epsilon = 1;
  result = equal(a,b,tau,epsilon);
  printf("And now a 0: %d\n",result);
  return 0;
}
