#include "bases.h"
#include <math.h>
double buildingvector_norm(matrix* m, int number_of_vector){
  double eredm;
  for (int i=0;i<m->rows;i++){
    eredm+=(extr_matrix_value(m, number_of_vector, i)*extr_matrix_value(m, number_of_vector, i));
  }
  return sqrt(eredm);
}


matrix * GramSchmidt (matrix* original){
  matrix * R = matrix_alloc(original->cols, original->cols);
  double* buffer_vector=malloc(sizeof(double)*original->rows);
  for (int i=0;i<original->cols;i++){
      for (int j=0;j<original->rows;j++){
          buffer_vector[j]=extr_matrix_value(original, i, j);//this is the ai "vector", the i. column of the A matrix
      }
    /*  if (i==0){
          for (int l=0;l<original->rows;l++){
              printf("buffer_vector %g \n", buffer_vector[l] );
          }
      }*/
      double Rii=buildingvector_norm(original,i);//i need this multiple times, so don't really want to calculate it separately
      add_matrix_value(R,i,i,Rii);
      for (int j=0;j<original->rows;j++){
          buffer_vector[j]/=Rii;//normalisation of the vector
          add_matrix_value(original, i, j, buffer_vector[j]);
      }
      for (int j=i+1;j<original->cols;j++){
          double Rij=0;
          for (int k=0;k<original->rows;k++){
              Rij+=extr_matrix_value(original, j, k)*buffer_vector[k];//this is the aj "vector", the j. column of the A matrix
        }
          add_matrix_value(R,j,i,Rij);
          add_matrix_value(R,i,j, 0);
    //      print_matrix(R);
          for (int k=0;k<original->rows;k++){
              double new_val=extr_matrix_value(original, j, k);
              new_val-=buffer_vector[k]*Rij;
              add_matrix_value(original, j, k, new_val);//this is the aj "vector", the j. column of the A matrix
          }
      }
  }
  return R;
}


vector*  qr_gs_solve(matrix* Q, matrix* R, vector* b){
  matrix* Qt=transpose(Q);
  vector* object=matrix_vector_appl(Qt, b);
  for (int i=b->n-1;i>-1;i--){
  //  printf("Rx(%d) %lg\n",i,object->data[i] );
    double buffer=0;
    for (int j=i+1;j<b->n;j++){
      buffer+=extr_matrix_value(R, j, i)*object->data[j];
  //    printf("Buffer: %lg\n", buffer );
    }
    object->data[i]=(object->data[i]-buffer)/extr_matrix_value(R, i,i);
//    printf("Rx(%d) now:%lg\n",i,object->data[i] );
  }
  matrix_free(Qt);
  return object;
}

int main() {
  matrix * m =random_matrix(4,4);
  matrix* A=copy_mtrx(m);
  matrix * e =matrix_alloc(4,4);
  printf("%i, %i\n", m->rows, m->cols);
  printf("original matrix\n");
  print_matrix(m);
  printf("Again\n" );
  matrix * mtr=matrix_multiplication(e,m);
  print_matrix(mtr);
  matrix * R =GramSchmidt(m);
  matrix * Qt=transpose(m);
  printf("...and Again\n" );
  mtr=matrix_multiplication(m, R);
  print_matrix(mtr);
  printf("matrix Q: \n \n" );
  print_matrix(m);

  printf("matrix R: \n \n" );
  print_matrix(R);
  printf("matrix Q transpose: \n \n" );
  print_matrix(Qt);
  matrix * egys=matrix_multiplication(Qt, m);
  printf("supposed to be eigenmatrix:\n" );
  print_matrix(egys);
  printf("Close enough\n" );
  printf("Now let's generate a random vector too:\n" );
  vector* b=vector_alloc(A->cols);
  for (int i=0;i<b->n;i++){
    b->data[i]=(double)rand()/(double)RAND_MAX;
    printf("%lg \t", b->data[i] );
  }
  printf("\n");
  vector* x=qr_gs_solve(m, R, b);
  printf("Now let's see the solution vector:\n" );
  for (int i=0;i<b->n;i++){
    printf("%lg \t", x->data[i]);
  }
  printf("\n");
  printf("and finally Ax = (supposed to be =b):\n" );
  vector* sol=matrix_vector_appl(A, x);
  for (int i=0;i<b->n;i++){
    printf("%lg \t", sol->data[i]);
  }
  printf("\n" );


    vector_free(x);
    vector_free(sol);
    vector_free(b);
    matrix_free(A);
    matrix_free(m);
    matrix_free(R);
    matrix_free(Qt);
    matrix_free(egys);
    matrix_free(e);
    matrix_free(mtr);


  printf(" Note that I never used the matrix's square-ness in the solution, and indeed, it would work on tall matrices too: \n" );
  matrix * m2 =random_matrix(6,4);
  matrix * e2 =matrix_alloc(6,6);
  printf("%i, %i\n", m2->rows, m2->cols);
  printf("original matrix\n");
  print_matrix(m2);
  printf("Again\n" );
  matrix * mtr2=matrix_multiplication(e2,m2);
  print_matrix(mtr2);
  matrix * R2 =GramSchmidt(m2);
  matrix * Qt2=transpose(m2);
  printf("...and Again\n" );
  mtr2=matrix_multiplication(m2, R2);
  print_matrix(mtr2);
  printf("matrix Q: \n \n" );
  print_matrix(m2);

  printf("matrix R: \n \n" );
  print_matrix(R2);
  printf("matrix Q transpose: \n \n" );
  print_matrix(Qt2);
  matrix * egys2=matrix_multiplication(Qt2, m2);
  printf("supposed to be eigenmatrix:\n" );
  print_matrix(egys2);
  printf("Close enough\n" );
  matrix_free(m2);
  matrix_free(R2);
  matrix_free(Qt2);
  matrix_free(egys2);
  matrix_free(e2);
  matrix_free(mtr2);


  return 0;
}
