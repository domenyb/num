#include "airy.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<math.h>


void vector_print(const char* s, gsl_vector* v){
  printf("%s",s);
  for(int i=0;i<v->size;i++){
    printf("%g",gsl_vector_get(v,i));
    printf("\t");
  }
  printf("\n");
}

void matrix_print(const char* s, gsl_matrix* A){
  printf("%s",s);
  for(int i=0;i<A->size1;i++){
    for(int j=0;j<A->size2;j++){
      printf("%g",gsl_matrix_get(A,i,j));
      printf("\t");
    }
  printf("\n");
  }
}

int main(){
	//int n, m;
	int size1, size2;

  scanf("%i",&size1);
  scanf("%i",&size2);
  printf("%i\t",size1);
  printf("%i\n",size2);

  gsl_matrix *A=gsl_matrix_alloc(size1,size2);
  gsl_matrix *B=gsl_matrix_alloc(size1,size2);
  gsl_matrix_fscanf(stdin,A);
  matrix_print("A: \n",A);

  gsl_vector *b=gsl_vector_alloc(size2);
  gsl_vector *x=gsl_vector_alloc(size2);

	gsl_vector_fscanf(stdin,b);
  vector_print("b: \n",b);
  gsl_matrix_memcpy(B,A);
  gsl_linalg_HH_solve(B,b,x);
  vector_print("x :\n",x);
  gsl_blas_dgemv(CblasNoTrans,1,A,x,0,b);
  vector_print("Ax : \n",b);
  printf("\n");

	airy_func();



gsl_matrix_free(A);
gsl_matrix_free(B);
gsl_vector_free(b);
gsl_vector_free(x);

return 0;
}
