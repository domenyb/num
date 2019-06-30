#include "sweep.h"
//#include "bases.h"
#include <math.h>


int Jacobi_cycles (matrix* A, vector* eigenvls, matrix* eigenvtrs){
  int sweep_nr=0;//from my cyclic jacobi excersise
  int more=1;
  for (int i=0;i<A->cols;i++){
    eigenvls->data[i]=extr_matrix_value(A,i,i);
  }
//  printf(" \n now the vector:\n" );
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
  matrix* A=random_simm_matrix(size);//this is a symmetric matrix
  matrix* D=random_simm_matrix_poz(size);//This is N, but this object will be D later - a positive definite symmetrical matrix
  printf("Original matrix A: \n \n" );
  print_matrix(A);
  printf("Original matrix N: \n \n" );
  print_matrix(D);
  matrix* N=copy_mtrx(D);//this is N where it keeps the value
  matrix* V=matrix_alloc(size,size);//my allocation function allocates an identity matrix by default
  matrix* E=matrix_alloc(size,size);//my allocation function allocates an identity matrix by default

  vector* not_needed=vector_alloc(size);//this is the eigenvalues of the N mtx, not needed in this excersise
  int sweep=Jacobi_cycles(D, not_needed, V);//decompose N to get V and D
  printf("Changed matrix D: \n \n" );
  print_matrix(D);

  matrix* Vt=transpose(V);
  matrix* sqrtD=sqrt_diag_mtx(D);//create sqrt D and the inverse of it
  printf("Changed matrix sqrt(D): \n \n" );
  print_matrix(sqrtD);
  matrix* sqrtDi=invert_diag_mtx(sqrtD);
  printf("Changed matrix sqrt(D)^(-1): \n \n" );
  print_matrix(sqrtDi);
  //now let's build B! for the multiplication part, I only need 2 buffer matrices to alternate in the in- and output
  matrix* buff1=matrix_alloc(size, size);
  matrix* buff2= matrix_alloc(size, size);
  matrix* B=matrix_alloc(size, size);
  matrix_multiplication_2(sqrtDi, V, buff1);//V * sqrt D^-1
  matrix_multiplication_2(buff1, A, buff2); // A * V * sqrt D^-1
  matrix_multiplication_2(buff2, Vt, buff1);// V^t * A * V * sqrt D^-1
  matrix_multiplication_2(buff1, sqrtDi, B);//B is the matrix we need
  printf("Now introduce the B matrix\n" );
  print_matrix(B);
  matrix_free(buff1);//i don't need the buffer matrices anymore
  matrix_free(buff2);
  printf("This is the matrix that needs the eigenvalues and vectors solved\n" );
  matrix* C=copy_mtrx(B);//save it for later control of the results
  vector* v=vector_alloc(size);//eigenvalues
  sweep=Jacobi_cycles(B, v, E);//E is the "V matrix of B"
  printf("The eigenvalues from smallest :\n" );
  for (int i=0;i<size;i++){
    printf("%g\t", v->data[i]);}//print the eigenvalues
  printf("\n" );
  matrix* Et=transpose(E);
  matrix* AV=matrix_multiplication(E,C);
  matrix* result=matrix_multiplication(AV,Et);// create the diagonal to verify the results
  printf("\n The diagonal matrix containing the above eigenvalues:\n \n" );
  print_matrix(result);
  printf("This is close enough to diagonal\n" );
  printf("Now we've found the eigenvalues (and vectors) of matrix B, to get the original eigenvectors: \n I should transform them back using the y= sqrt(D)*V^T x equation, where y, D, and V^T are known. \n" );
  //V is the eigenvtrs

  //check whether the vales we got are the correct ones
  vector * eigen1=eigen_vector(E, 0);
  printf("First eigenvector:\n" );
  for (int i=0;i<eigen1->n;i++){
    printf("%g \t" , eigen1->data[i] );
  }
  printf("\n" );
  vector * eig_multip=matrix_vector_appl(C, eigen1);
  printf("The matrix applied to the first eigenvector:\n" );
  for (int i=0;i<eigen1->n;i++){
    printf("%g \t" , eig_multip->data[i] );
  }
  printf("\n" );
  printf("The first eigenvector multiplied by the eigenvalue:\n" );
  for (int i=0;i<eigen1->n;i++){
    printf("%g \t" , eigen1->data[i]*v->data[0]);
  }
  printf("\n The two are identical, so we really got the good values. Transforming all the y vectors back to x: \n first use sqrt(D)^-1, then V to make the two Transforming matrices back to identity matrices\n" );

  for (int i=0;i<A->rows;i++){
    vector* eigv_tr1=eigen_vector(E, i);
    vector* eigv_tr2=matrix_vector_appl(sqrtDi, eigv_tr1);
    vector* eigv=matrix_vector_appl(V, eigv_tr2);
    printf("vector %d: \n",i );
    for (int j=0;j<eigv->n;j++){
      printf("%g\t", eigv->data[j]  );
    }
    printf("\n \n" );
    vector_free(eigv_tr2);
    vector_free(eigv_tr1);
    vector_free(eigv);//free the stuff up created in the loop
  }
  printf("This concludes the exam excersise. The required analytical proof is in a separate pdf file (out.pdf)\n" );
  //free up all the space allocated - double checked, I do have all of them correctly



  matrix_free(B);
  matrix_free(A);
  matrix_free(D);
  matrix_free(E);
  matrix_free(Et);
  matrix_free(AV);
  matrix_free(V);
  matrix_free(Vt);
  matrix_free(result);
  vector_free(not_needed);
  vector_free(v);
  matrix_free(N);
  matrix_free(sqrtDi);
  matrix_free(sqrtD);
  matrix_free(C);
  vector_free(eigen1);
  vector_free(eig_multip);
  return 0;
}
