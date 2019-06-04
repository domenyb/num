#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"funcs.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>


int main(int argc, char** argv){
  //from the gsl documentation
  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter = 0;
  int status;
  double size;

  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 0.0);
  gsl_vector_set (x, 1, 0.0);

  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  minex_func.n = 2;
  minex_func.f = rosen_func;
  minex_func.params = NULL;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      printf ("%5ld %10.3e %10.3e The Rosen. func:%7.3f size:%.3f\n",
              iter,
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

//second part
  double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
  double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
  double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
  int n = sizeof(t)/sizeof(t[0]);
  int dim=3;
  struct experimental_data data;
  data.n=n;
  data.t=t;
  data.y=y;
  data.e=e;
  gsl_multimin_function F;
  F.f = function_to_minimize;
  F.n = dim;
  F.params=(void*)&data;

  gsl_multimin_fminimizer *M;
  #define TYPE gsl_multimin_fminimizer_nmsimplex2
  M = gsl_multimin_fminimizer_alloc(TYPE,dim);
  gsl_vector* start=gsl_vector_alloc(dim);
  gsl_vector* step=gsl_vector_alloc(dim);
  gsl_vector_set(start,0,6);
  gsl_vector_set(start,1,4);
  gsl_vector_set(start,2,1);
  gsl_vector_set(step,0,0.2);
  gsl_vector_set(step,1,0.2);
  gsl_vector_set(step,2,0.1);

  gsl_multimin_fminimizer_set(M,&F,start,step);

  int iter2=0,status2=GSL_CONTINUE;
  double size2;
  while (status2 == GSL_CONTINUE && iter2 < 100){
    iter2++;
    status2 = gsl_multimin_fminimizer_iterate(M);
    if (status2) break;

    size2 = gsl_multimin_fminimizer_size (M);
    status2 = gsl_multimin_test_size (size2, 1e-2);


      printf ("iteration:%5d \t A=%g \t  T=%g \t  B=%g \t  function:%g \t  size:%g\n", iter2,gsl_vector_get (M->x, 0), gsl_vector_get (M->x, 1), gsl_vector_get (M->x, 2), M->fval, size2);
    }


  double A=gsl_vector_get(M->x,0);
  double Tc=gsl_vector_get(M->x,1);
  double B=gsl_vector_get(M->x,2);
  for(int i=0;i<n;i++) {
    fprintf(stderr,"%g %g %g %g\n",*(t+i),*(y+i),*(e+i),decayfunction(t[i],A,Tc,B)
  );}


  printf("A=%g \t T=%g \t B=%g\n",A,Tc,B);


  return 0;
}
