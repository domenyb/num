#include <stdio.h>
#include<stdlib.h>
typedef struct {int n; double * data;} vector;
typedef struct {int rows,cols; double * data;} matrix;



vector* vector_alloc(int length){
  vector* object=(vector*) malloc(sizeof(vector));
  object->n=length;
  object->data=malloc(length*sizeof(double));
  for (int i=0;i<length;i++){
    object->data[i]=0;
  }
  return object;
}

void vector_free(vector* v){
  free(v->data);
  free(v);
}


void add_matrix_value(matrix* m, int col, int row, double value) {
  int place=m->cols*row+col;//i wouldn't need to calculate this separately, but this is "cleaner" for practice code
  m->data[place]=value;
  }

void matrix_reset(matrix* m){
  for (int i=0;i<m->cols;i++){
    for (int j=0;j<m->rows;j++){
      if (i==j) add_matrix_value(m, i,j,1);
      else add_matrix_value(m, i,j,0);
    }
  }
}

matrix * matrix_alloc(int rows, int cols){
  //also builds it to be an identity matrix (or sort of...)
  matrix * object =(matrix *) malloc(sizeof(matrix));
  object->rows=rows;
  object->cols=cols;
  object->data=malloc(rows*cols*sizeof(double));
  for (int i=0;i<rows;i++){
    for (int j=0;j<cols;j++){
      if (i==j){
        add_matrix_value(object, j, i, 1);
      }
      else {
        add_matrix_value(object, j, i, 0);
      }
    }
  }
  return object;
}

void matrix_free(matrix * m){
  free(m->data);
  free(m);
}

double extr_matrix_value(matrix* m, int col, int row) {
  int place=m->cols*row+col;//i wouldn't need to calculate this separately, but this is "cleaner" for practice code
  return m->data[place];
  }

void print_matrix(matrix * m){
//  printf("Megkaptam a matrixot: %i, %i\n", m->rows, m->cols);
  for (int i=0;i<m->rows;i++){
//    printf("Row_leap\n" );
    for (int j=0;j<m->cols;j++){
      printf("%g\t", extr_matrix_value(m, j, i));
    }
    printf("\n" );
  }
}

matrix * read_matrix(FILE* stored){
  int rows, cols;
  double buff;
  fscanf(stored, "%d", &rows);
  fscanf(stored, "%d", &cols);
  matrix * object=matrix_alloc(rows, cols);
  for (int i=0;i<rows;i++){
//    printf("Row_leap\n" );
    for (int j=0;j<cols;j++){
      fscanf(stored, "%lg", &buff);
      add_matrix_value(object, j, i, buff);
    }
  }
    return object;
}

matrix * matrix_multiplication(matrix * m1, matrix * m2){
  //I realise I should at first check whether the matrices are multiplicable, but this function will likely never see incompatible matrices, since it's purpose is to verify the others
  matrix * object=matrix_alloc(m1->rows, m2->cols);
  for (int k=0;k<m2->cols;k++){
//    printf("k=%d\n", k);
    for (int i=0;i<m1->rows;i++){
      double buff=0;
//      printf("i=%d\n", i);
      for (int j=0;j<m1->cols;j++){
  //      printf("j=%d\n", j);
        buff+=extr_matrix_value(m1, j, i)*extr_matrix_value(m2, k, j);
      }
      add_matrix_value(object, k,i, buff);
    }
  }
  return object;
}
matrix* transpose(matrix* m){
  matrix * object=matrix_alloc(m->cols, m->rows);
  for (int i=0;i<m->cols;i++){
    for (int j=0;j<m->rows;j++){
      add_matrix_value(object, j, i, extr_matrix_value(m, i, j));
    }
  }
  return object;
}


matrix * random_matrix(int rows, int cols){
  matrix * object=matrix_alloc(rows, cols);
  for (int i=0;i<rows;i++){
//    printf("Row_leap\n" );
    for (int j=0;j<cols;j++){
      add_matrix_value(object, j, i, (double)rand()/(double)RAND_MAX);
    }
  }
    return object;
}
matrix* copy_mtrx(matrix* m){

    matrix * object=matrix_alloc(m->rows, m->cols);
    for (int i=0;i<m->cols;i++){
      for (int j=0;j<m->rows;j++){
        add_matrix_value(object, i, j, extr_matrix_value(m, i, j));
      }
    }
    return object;
  }

matrix * random_simm_matrix(int size){
  matrix* object=matrix_alloc(size, size);
  for (int i=0;i<size;i++){
    for (int j=i;j<size;j++){
      double buffer=(double)rand()/(double)RAND_MAX;
      add_matrix_value(object, j, i, buffer);
      add_matrix_value(object, i, j, buffer);
    }
  }
  return object;
}

vector* matrix_vector_appl(matrix* m, vector* v){
  vector* object=vector_alloc(v->n);
  for (int i=0;i<v->n;i++){
    double buffer=0;
    for (int j=0;j<v->n;j++){
      buffer+=extr_matrix_value(m, j, i)*v->data[j];
    }
    object->data[i]=buffer;
  }
  return object;
}
