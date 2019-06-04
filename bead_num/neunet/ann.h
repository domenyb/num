#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multimin.h>
#include<gsl/gsl_blas.h>
#include<math.h>
//after the example codes - now I want to use the gsl vector functions, since the bases.h doesn't host quasinewtonian - and that excersise isn't done by me at the time of making this

typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;
ann* ann_alloc(int n, double(*f)(double)){
  ann* object = malloc(sizeof(ann));
  object->n = n;
  object->f = f;
  object->data= gsl_vector_alloc(3*n);
  return object;
}

void ann_free(ann* network){
  gsl_vector_free(network->data);
  free(network);
}

double ann_feed_forward(ann* network, double x){
  double s=0;
  double a, b, w, y;
  for (int i=0; i<network->n;i++) {
    a=gsl_vector_get(network->data,3*i);
    b=gsl_vector_get(network->data,3*i+1);
    w=gsl_vector_get(network->data,3*i+2);
    y=network->f((x-a)/b)*w/b;
    s+=y;
  }
  return s;
  }

double ann_feed_forwardB(ann* network, double x, double (*derivatives)(double)){//for the derivative - basically I throw the same things through the neunet, but now I don't use the neuron's own function, but the derivative
  double s=0;
  double a, b, w, y;
  for (int i=0;i<network->n;i++){
    a=gsl_vector_get(network->data,3*i);
    b=gsl_vector_get(network->data,3*i+1);
    w=gsl_vector_get(network->data,3*i+2);
    y=derivatives((x-a)/b)*w/b;//I could load either the derivative or the antiderivate of my function (both are ugly)
    s+=y;
  }
  return s;
}




void ann_train(ann* network, gsl_vector* vx, gsl_vector* vy){
  double delta(const gsl_vector* p, void* params){//apparently this is needed for the gsl multimin....
      gsl_vector_memcpy(network->data,p);
      double result=0;
  		for(int i=0; i<vx->size; i++){
  			double x=gsl_vector_get(vx,i);
  			double y=gsl_vector_get(vy,i);
  			double f=ann_feed_forward(network,x);
  			result+=fabs(f-y);
  		}
  		return result/vx->size;
  	}
  //minimizing after the minimization function of the gsl
  gsl_vector* p=gsl_vector_alloc(network->data->size);
  gsl_vector_memcpy(p,  network->data);
  gsl_multimin_function minex;
  int iter=0, status; double size;
  minex.f=delta;
  minex.n=network->data->size;
  minex.params=NULL;
  gsl_vector *ss = gsl_vector_alloc(minex.n);
  gsl_vector_set_all(ss,0.5);
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,minex.n);
  gsl_multimin_fminimizer_set (s, &minex, p, ss);
  do{
  		iter++;
  		status = gsl_multimin_fminimizer_iterate(s);
  		if(status){
  			fprintf(stderr,"done? \n");
  			break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);
  	}while( status == GSL_CONTINUE && iter < 10000);
  gsl_vector_memcpy(network->data,s->x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(p);
}
