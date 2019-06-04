#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"bases.h"



void testfunc(double t, vector* y, vector*  dydx){
//I mean... this is an easy thing to think about - an evergreen example
  dydx->data[0]=y->data[1];
  double y_1 = y->data[0];
  dydx->data[1]=1-y_1+0.01*y_1*y_1;
}



//#include"odefunctions.h"
void rkstep12(	double t,                                  /* the current value of the variable */
	double h,                                  /* the step to be taken */
	vector* yt,                                /* the current value y(t) of the sought function */
	void f(double t, vector* y, vector* dydt), /* the right-hand-side, dydt = f(t,y) */
	vector* yth,                               /* output: y(t+h) */
	vector* err                                /* output: error estimate dy */
){//after the one in the lecture notes
//  vector* pszeu=vector_alloc(yt->n);
  f(t, yt, err);

  for (int i=0;i<err->n;i++){
    err->data[i]*=h/2;
    yt->data[i]+=err->data[i];
  }
f(t+h/2, yt, yth);

  for (int i=0;i<err->n;i++){
    yt->data[i]-=err->data[i];
    err->data[i]+=yth->data[i]*(-h/2.0);
    yth->data[i]*=h;
    yth->data[i]+=yt->data[i];
  }
}

void driver(
	double* t,                             /* the current value of the variable */
	double b,                              /* the end-point of the integration */
	double* h,                             /* the current step-size */
	vector*yt,                             /* the current y(t) */
	double acc,                            /* absolute accuracy goal */
	double eps,                            /* relative accuracy goal */
	void stepper(                          /* the stepper function to be used */
		double t, double h, vector*yt,
		void f(double t,vector*y,vector*dydt),
		vector*yth, vector*err
		),
	void f(double t,vector*y,vector*dydt) /* right-hand-side */
){//following the lecture notes
  int sweep=0;
  int dim=yt->n;
  vector* ycurr=vector_alloc(dim);
  vector* err=vector_alloc(dim);
  double starting_point=*t;
  double norm1;
  double norm2;
  double tolerance;
  while(fabs(*t - b) > 1e-12 && sweep < 1000000){
      if(fabs(*t+*h)>/*< no, other way around*/ fabs(b)){
        *h = b - *t;
      }
      stepper(*t, *h, yt, f, ycurr, err);
      norm1=vector_norm(err);
      norm2=vector_norm(yt);
      tolerance=(acc+eps*norm2)*sqrt(*h/(b-starting_point));
      if (tolerance>norm1){
        sweep++;

        for (int i=0;i<dim;i++){
          yt->data[i]=ycurr->data[i];
        }
        *t+=*h;
      }
      //update the error
    *h*=pow(tolerance/norm1, 0.25)*0.95;
  }

  vector_free(ycurr);
  vector_free(err);

}


void driver_B(
	double* t,                             /* the current value of the variable */
	double b,                              /* the end-point of the integration */
	double* h,                             /* the current step-size */
	vector*yt,                             /* the current y(t) */
	double acc,                            /* absolute accuracy goal */
	double eps,                            /* relative accuracy goal */
	void stepper(                          /* the stepper function to be used */
		double t, double h, vector*yt,
		void f(double t,vector*y,vector*dydt),
		vector*yth, vector*err
		),
	void f(double t,vector*y,vector*dydt) /* right-hand-side */
  , matrix* path, int* i_am_here //path: i, j, let i be the iteration nr, j the coordinate
){//following the lecture notes + modification
  int sweep=0;
  int dim=yt->n;
  vector* ycurr=vector_alloc(dim);
  vector* err=vector_alloc(dim);
  double starting_point=*t;
  double norm1;
  double norm2;
  double tolerance;
  while(fabs(*t - b) > 1e-12 && sweep < 1000000){
      if(fabs(*t+*h)>/*< no, other way around*/ fabs(b)){
        *h =b-*t;
      }
      stepper(*t, *h, yt, f, ycurr, err);
      norm1=vector_norm(err);
      norm2=vector_norm(yt);
      tolerance=(acc+eps*norm2)*sqrt(*h/(b-starting_point));
      if (tolerance>norm1){
        sweep++;
        for (int i=0;i<dim;i++){
          yt->data[i]=ycurr->data[i];

        }
//        printf(" %g \t %g  \t %g\n", *t, yt->data[0], yt->data[1] );
        add_matrix_value(path, 0, *i_am_here, *t);
        for (int i=0;i<yt->n;i++){
          add_matrix_value(path, i+1, *i_am_here, yt->data[i]);
        }
        *t+=*h;
        (*i_am_here)++;
      }
      //update the error
    *h*=pow(tolerance/norm1, 0.25)*0.95;
  }

  vector_free(ycurr);
  vector_free(err);
}
