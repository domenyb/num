#include "sweep.h"
//#include "bases.h"
#include <math.h>


int Jacobi_cycles (matrix* A, vector* eigenvls, matrix* eigenvtrs, int smallest,int num_rows){
  int sweep_nr=0;
  int more=1;
  for (int i=0;i<A->cols;i++){
    eigenvls->data[i]=extr_matrix_value(A,i,i);
  }
//  printf(" \n now the vector:\n" );
  for (int i=0;i<eigenvls->n;i++){
//    printf(" %g\t", eigenvls->data[i]);
  }
 // printf("\n" );
//  for (int i=0;i<num_rows;i++){
  while (more !=0){
    more=jacobi_sweep_one_row(A, eigenvls, eigenvtrs, num_rows, smallest);//pipa
    sweep_nr++;
//    printf("Sweep %d\n", sweep_nr );printf("%d : %g\n",i, eigenvls->data[i] );
  //  if (sweep_nr>100) return 0;
//}
}
  return sweep_nr;
}


int main(int argc, char **argv ) {
  int size=atoi(argv[1]);
  int numrows=atoi(argv[2]);
  int smallest=1;
  if (argc>3){
    if (atoi(argv[3])==0) smallest=0;
  }

  matrix* A=random_simm_matrix(size);
  printf("Original matrix: \n \n" );
  print_matrix(A);
  matrix* B=copy_mtrx(A);
  matrix* E=matrix_alloc(size,size);
  printf("Eigenmatrix: \n \n" );
  print_matrix(E);
  vector* v=vector_alloc(size);
  int sweep=Jacobi_cycles(A, v, E, smallest, numrows);
  printf("Changed matrix: \n \n" );
  print_matrix(A);
//  printf("# of sweeps: %d \n", sweep);
  if (smallest!=0){
  printf("The first %d eigenvalues from smallest :\n", numrows );}
  else {printf("The first %d eigenvalues from largest :\n", numrows );}
  for (int i=0;i<numrows;i++){
    printf("%g\t", v->data[i]);
  }
  printf("\n" );
  matrix* Vt=transpose(E);
  matrix* AV=matrix_multiplication(E,B);
  matrix* result=matrix_multiplication(AV,Vt);
/*  printf("\n The diagonal matrix containing the above eigenvalues:\n \n" );
  print_matrix(result);*/
  printf("Close enough to diagonal\n" );
  matrix_free(B);
  matrix_free(A);
  matrix_free(E);
  matrix_free(Vt);
  matrix_free(result);
  vector_free(v);

  return 0;
}
