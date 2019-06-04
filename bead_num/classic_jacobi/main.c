#include "sweep.h"
//#include "bases.h"
#include <math.h>




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
  //if (smallest!=0){
  printf("The eigenvalues from smallest :\n" );
  //else {printf("The first %d eigenvalues from largest :\n", numrows );}
  for (int i=0;i<size;i++){
    printf("%g\t", v->data[i]);}
  //}
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
