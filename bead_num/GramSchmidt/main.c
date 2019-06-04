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

matrix* invert_upper_triangular(matrix* R){
  matrix* object=matrix_alloc(R->rows, R->cols);
  for (int i=0;i<R->rows;i++){
    add_matrix_value(object, i, i, 1/extr_matrix_value(R, i, i));
  }
  for (int i=0;i<R->rows;i++){
    for (int j=i+1;j<R->rows;j++){
      double buffer=0;
      for (int k=0;k<j;k++){
        buffer-=extr_matrix_value(object, k, i)*extr_matrix_value(R, j, k);
//        printf("%g*%g=%g\n",extr_matrix_value(object, k, i), extr_matrix_value(R, j, k), buffer);
      }
      buffer=buffer/extr_matrix_value(R, j,j);
  //    printf("%g\n", buffer);
      add_matrix_value(object, j, i, buffer);
//      print_matrix(object);
//      printf("\n");
    }
  }
  return object;

}


matrix*  qr_gs_inverse( matrix* Q, matrix* R){
    matrix* Qt=transpose(Q);
    matrix* Rt=invert_upper_triangular(R);
    matrix * object=matrix_multiplication(Rt, Qt);
    matrix_free(Qt);
    matrix_free(Rt);
    return object;
}

int main() {
  //generate random matrix
  matrix * m =random_matrix(10,10);
  matrix* A=copy_mtrx(m);
  printf("original matrix\n");
  print_matrix(m);/*
  printf("Again\n" );
  matrix * mtr=matrix_multiplication(e,m);
  print_matrix(mtr);*/
  matrix * R =GramSchmidt(m);
  matrix * Qt=transpose(m);
//  printf("...and Again\n" );
//  mtr=matrix_multiplication(m, R);
//  print_matrix(mtr);
  matrix* B=qr_gs_inverse(m, R);
  printf("matrix Q: \n \n" );
  print_matrix(m);
  printf("matrix R: \n \n" );
  print_matrix(R);/*
  printf("matrix Rt: \n \n" );
  print_matrix(Rt);*/

  printf("matrix B: \n \n" );
  print_matrix(B);
  matrix * egys=matrix_multiplication(A, B);
  printf("\n \n supposed to be eigenmatrix:\n" );
  print_matrix(egys);
  matrix_free(m);
  matrix_free(R);
  matrix_free(Qt);
  matrix_free(egys);
  matrix_free(B);
//  matrix_free(e);
//  matrix_free(mtr);
  printf("Close enough\n" );
  return 0;
}
