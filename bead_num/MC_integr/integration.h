#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bases.h"
#include <assert.h>
#define rnd (double)rand()/(double)RAND_MAX

void randomx(vector* a, vector* b, vector* x){
  for (int i=0;i<a->n; i++){
    double buff=rnd;//this is a number between 0 and 1, a+rnd(b-a) is sure to be between a and b, which is what we want
    x->data[i]=a->data[i]+buff*(b->data[i]-a->data[i]);
  }
}


void plainmc(double f(vector* x),vector* a, vector * b, int N, double* result, double* error){//after the text of the excersise
//  randomx = function(a,b) [a[i]+Math.random()*(b[i]-a[i]) for (i in a)]; I actually took this out, don't want it to be nested
  double volume=1;
  for(int i=0;i<a->n;i++)volume*=b->data[i]-a->data[i];
  double sum=0,sum2=0;
  vector* x=vector_alloc(a->n);
  for(int i=0;i<N;i++){
    randomx(a,b,x);
    double fx=f(x);
    sum += fx;
    sum2 += fx*fx;
  }
  double mean = sum/N;                             // <f_i>
  double sigma = sqrt(sum2/N-mean*mean);    // sigma² = <(f_i)²> - <f_i>²
  double SIGMA = sigma/sqrt(N);               // SIGMA² = <Q²> - <Q>²
  *result= mean*volume;
  *error=SIGMA*volume;
  vector_free(x);
}
