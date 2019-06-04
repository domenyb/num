#include "sweep.h"
//#include "bases.h"
#include <math.h>


int Jacobi_cycles (matrix* A, vector* eigenvls, matrix* eigenvtrs){
  int sweep_nr=0;
  int more=1;
  for (int i=0;i<A->cols;i++){
    eigenvls->data[i]=extr_matrix_value(A,i,i);
  }
  printf(" \n now the vector:\n" );
  for (int i=0;i<eigenvls->n;i++){
    printf(" %g\t", eigenvls->data[i]);
  }
  printf("\n" );
  while (more !=0){
    more=jacobi_sweep(A, eigenvls, eigenvtrs);//pipa
    sweep_nr++;
    printf("Sweep %d\n", sweep_nr );
  }
  return sweep_nr;
}


int main(int argc, char **argv ) {
  int size=atoi(argv[1]);
  matrix* A=random_simm_matrix(size);
  printf("Original matrix: \n \n" );
  print_matrix(A);
  matrix* B=copy_mtrx(A);
  matrix* E=matrix_alloc(size,size);
  printf("Eigenmatrix: \n \n" );
  print_matrix(E);
  vector* v=vector_alloc(size);
  int sweep=Jacobi_cycles(A, v, E);
  printf("Changed matrix: \n \n" );
  print_matrix(A);
  printf("# of sweeps: %d \n", sweep);
  printf("The eigenvalues:\n" );
  for (int i=0;i<v->n;i++){
    printf("%g\t", v->data[i]);
  }
  printf("\n" );
  matrix* Vt=transpose(E);
  matrix* AV=matrix_multiplication(E,B);
  matrix* result=matrix_multiplication(AV,Vt);
  printf("\n The diagonal matrix containing the above eigenvalues:\n \n" );
  print_matrix(result);
  printf("Close enough to diagonal\n" );
  matrix_free(B);
  matrix_free(A);
  matrix_free(E);
  matrix_free(Vt);
  matrix_free(result);
  vector_free(v);

  return 0;
}
