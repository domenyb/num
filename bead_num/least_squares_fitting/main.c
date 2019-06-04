#include "root.h"

double calc_func(double funs(int i, double x),int funsnr, vector* params, vector* errors, int sgn, double x){//sgn should be +1/-1/0 depending on the values
  double buffer=0;
  for (int i=0;i<funsnr;i++){
    buffer+=funs(i,x)*(params->data[i]+sgn*errors->data[i]);
  }
  return buffer;
}


int main(int argc, char** argv){
  int freeprm=3;
  if (argc==1) {printf("Need argument: input file name\n" ); return 1;}
  if (argc>2) {printf("Need only one argument: program <input file name>\n" ); return 1;}
  FILE* f=fopen(argv[1], "r");
  int length;
  fscanf(f, "%d", &length);
  vector* x=vector_read(f, length);
  vector* y=vector_read(f, length);
  vector* dy=vector_read(f, length);
  fclose(f);
  matrix* covariance=matrix_alloc(freeprm, freeprm);
  vector* coeffs=least_squares_fit(funs, freeprm, x, y, dy, covariance );

  printf("Fitting done, parameters:\n" );
  for (int i=0;i<freeprm;i++){
    printf("%lg\t", coeffs->data[i]);
  }
  printf("\n" );
  printf("covariance matrix:\n" );
  print_matrix(covariance);
  //from the cov matrix, let's extract the error values of C vector
  vector* errors=vector_alloc(freeprm);
  for (int i=0;i<freeprm;i++){
    errors->data[i]=sqrt(extr_matrix_value(covariance, i,i));
  }

  printf(" giving out the fitting values\n" );

  FILE* orig_data=fopen("out1.dat", "w");
  FILE* fitted=fopen("out2.dat", "w");
  FILE* fit_perr=fopen("out3.dat", "w");
  FILE* fit_merr=fopen("out4.dat", "w");

  for (int i=0;i<length;i++){
    fprintf(orig_data, "%lg \t %lg \t %lg \n",x->data[i],y->data[i],dy->data[i] );
  }

  for (double i=0;i<10;i+=0.01){
    fprintf(fitted, "%lg \t %lg \n", i, calc_func(funs, freeprm, coeffs, errors, 0, i) );
    fprintf(fit_perr, "%lg \t %lg \n", i, calc_func(funs, freeprm, coeffs, errors, +1, i) );
    fprintf(fit_merr, "%lg \t %lg \n", i, calc_func(funs, freeprm, coeffs, errors, -1, i) );
  }
  vector_free(coeffs);
  vector_free(x);
  vector_free(y);
  vector_free(dy);
  vector_free(errors);
  matrix_free(covariance);
  fclose(orig_data);
  fclose(fitted);
  fclose(fit_perr);
  fclose(fit_merr);

  return 0;
}
