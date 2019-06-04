#include "bases.h"

void func1(const vector* x, vector* fx, matrix* jacobian){
  double A=10000;
  fx->data[0]=A*x->data[0]*x->data[1]-1;//  A*x*y = 1 ,
  fx->data[1]=exp(-1*x->data[0])+exp(-1*x->data[1])-1-1/A;//exp(-x) + exp(-y) = 1 + 1/A,
  double j11 = A * x->data[1] ; // df1/dx
  double j12 = A * x->data[0]; // df1/dy
  double j21 = -exp(-1*x->data[0]); // df2/dx
  double j22 = -exp(-1*x->data[1]); // df2/dy
  add_matrix_value(jacobian, 0,0,j11);
  add_matrix_value(jacobian, 0,1,j12);
  add_matrix_value(jacobian, 1,0,j21);
  add_matrix_value(jacobian, 1,1,j22);


  //we don't need to return anything


}


void func2(const vector* x, vector* fx, matrix* jacobian){
  //the name of these func-s are ugly, anyways solution should be (a,a^2) where in our case a=1, so (1,1)
  double x1=x->data[0];
  double x2=x->data[1];
  //f(x,y) = (1-x)^2+100(y-x^2)^2 is the function, gradient: df/dx=-2*(1-x)+100*2*(y-x^2)*(-2x), and df/dy=200(y-x^2)
  fx->data[0]=-2*(1-x1)+100*2*(x2-x1*x1)*-2*x1;
  fx->data[1]=200*(x2-x1*x1);
  double j11=2+400*x2+1200*x1*x1;
  double j12=-400*x1;
  double j21=-400*x1;//it would be strange if these two wouldn't match
  double j22=200;
  add_matrix_value(jacobian, 0,0,j11);
  add_matrix_value(jacobian, 0,1,j12);
  add_matrix_value(jacobian, 1,0,j21);
  add_matrix_value(jacobian, 1,1,j22);
}


void func3(const vector* x, vector* fx, matrix* jacobian){/*It has one local maximum at x = − 0.270845 { x=-0.270845} { x=-0.270845} and y = − 0.923039 { y=-0.923039} { y=-0.923039} where f ( x , y ) = 181.617 { f(x,y)=181.617} { f(x,y)=181.617}, and four identical local minima:

    f ( 3.0 , 2.0 ) = 0.0 , { f(3.0,2.0)=0.0, } f(3.0,2.0)=0.0,
    f ( − 2.805118 , 3.131312 ) = 0.0 , { f(-2.805118,3.131312)=0.0, } f(-2.805118,3.131312)=0.0,
    f ( − 3.779310 , − 3.283186 ) = 0.0 , { f(-3.779310,-3.283186)=0.0, } f(-3.779310,-3.283186)=0.0,
    f ( 3.584428 , − 1.848126 ) = 0.0. { f(3.584428,-1.848126)=0.0. } f(3.584428,-1.848126)=0.0.
    from wikipedia. We should find exactly one of these
  */
  double x1=x->data[0];
  double x2=x->data[1];
  //f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2, gradient: df/dx=2*2*(x^2+y-11)*x+2*(x+y^2-7) df/dy=2(x^2+y-11)+2*2*(x+y^2-7)*y
  fx->data[0]=4*(x1*x1 + x2 - 11)*x1 + 2*(x1 + x2*x2 - 7);
  fx->data[1]=2*(x1*x1 + x2 - 11) + 4*(x1 + x2*x2 - 7)*x2;
  double j11= 12*x1*x1 + 4*x2-42; //printf("j11 = %g \n", j11 );
  double j12 = 4*(x1+x2);
  double j21 = 4*(x1+x2);
  double j22 = 4*x1+12*x2*x2-26;
  add_matrix_value(jacobian, 0,0,j11);
  add_matrix_value(jacobian, 0,1,j12);
  add_matrix_value(jacobian, 1,0,j21);
  add_matrix_value(jacobian, 1,1,j22);

}


void func1_noJ(const vector* x, vector* fx){
  double A=10000;
  fx->data[0]=A*x->data[0]*x->data[1]-1;//  A*x*y = 1 ,
  fx->data[1]=exp(-1*x->data[0])+exp(-1*x->data[1])-1-1/A;//exp(-x) + exp(-y) = 1 + 1/A,

  //we don't need to return anything


}


void func2_noJ(const vector* x, vector* fx){
  //the name of these func-s are ugly, anyways solution should be (a,a^2) where in our case a=1, so (1,1)
  double x1=x->data[0];
  double x2=x->data[1];
  //f(x,y) = (1-x)^2+100(y-x^2)^2 is the function, gradient: df/dx=-2*(1-x)+100*2*(y-x^2)*(-2x), and df/dy=200(y-x^2)
  fx->data[0]=-2*(1-x1)+100*2*(x2-x1*x1)*-2*x1;
  fx->data[1]=200*(x2-x1*x1);
}


void func3_noJ(const vector* x, vector* fx){/*It has one local maximum at x = − 0.270845 { x=-0.270845} { x=-0.270845} and y = − 0.923039 { y=-0.923039} { y=-0.923039} where f ( x , y ) = 181.617 { f(x,y)=181.617} { f(x,y)=181.617}, and four identical local minima:

    f ( 3.0 , 2.0 ) = 0.0 , { f(3.0,2.0)=0.0, } f(3.0,2.0)=0.0,
    f ( − 2.805118 , 3.131312 ) = 0.0 , { f(-2.805118,3.131312)=0.0, } f(-2.805118,3.131312)=0.0,
    f ( − 3.779310 , − 3.283186 ) = 0.0 , { f(-3.779310,-3.283186)=0.0, } f(-3.779310,-3.283186)=0.0,
    f ( 3.584428 , − 1.848126 ) = 0.0. { f(3.584428,-1.848126)=0.0. } f(3.584428,-1.848126)=0.0.
    from wikipedia. We should find exactly one of these
  */
  double x1=x->data[0];
  double x2=x->data[1];
  //f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2, gradient: df/dx=2*2*(x^2+y-11)*x+2*(x+y^2-7) df/dy=2(x^2+y-11)+2*2*(x+y^2-7)*y
  fx->data[0]=4*(x1*x1 + x2 - 11)*x1 + 2*(x1 + x2*x2 - 7);
  fx->data[1]=2*(x1*x1 + x2 - 11) + 4*(x1 + x2*x2 - 7)*x2;
}



double buildingvector_norm(matrix* m, int number_of_vector){
  double eredm;
  for (int i=0;i<m->rows;i++){
    eredm+=(extr_matrix_value(m, number_of_vector, i)*extr_matrix_value(m, number_of_vector, i));
  }
  return sqrt(eredm);
}



matrix*  qr_gs_inverse( matrix* Q, matrix* R){
    matrix* Qt=transpose(Q);
    matrix* Rt=invert_upper_triangular(R);
    matrix * object=matrix_multiplication(Rt, Qt);
    matrix_free(Qt);
    matrix_free(Rt);
    return object;
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
//    printf("Rx(%d) %lg\n",i,object->data[i] );
    double buffer=0;
    for (int j=i+1;j<b->n;j++){
      buffer+=extr_matrix_value(R, j, i)*object->data[j];
//      printf("Buffer: %lg\n", buffer );
    }
    object->data[i]=(object->data[i]-buffer)/extr_matrix_value(R, i,i);
//    printf("Rx(%d) now:%lg\n",i,object->data[i] );
  }
  matrix_free(Qt);
  return object;
}

double funs(int i, double x){
   switch(i){
   case 0: return log(x); break;
   case 1: return 1.0;   break;
   case 2: return x;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

vector* least_squares_fit(double funs(int i, double x), int free_par, vector* x_data, vector* y_data, vector* dy_data, matrix* covariance) {
  int n=x_data->n;

  matrix* A = matrix_alloc(n,free_par);
  vector* b = vector_alloc(n);
  matrix* I = matrix_alloc(free_par, free_par); //identity matrix generated by default

  for (int i=0;i<n;i++){
    b->data[i]=y_data->data[i]/dy_data->data[i];
    for(int k=0;k<free_par;k++){
      add_matrix_value(A, k, i, funs(k, x_data->data[i])/dy_data->data[i]);
    }
  }
  matrix* R=GramSchmidt(A);
  //now A is Q
  vector* sol=qr_gs_solve(A,R,b);//solution vector - the three coefficients
  //covariance matrix? cov=R inverse* R inverse^transpose
  //for that, I'll need the inverse of R, which we fortunately do have
  //R is an upper triangular matrix, and we have an "invert upper triangular matrix" function from back at the Gram Schmidt excersise
  matrix* R_inverse=invert_upper_triangular(R);
  matrix* RiT=transpose(R_inverse);
  //here comes the problem, that all my funcs are not void, looking back that would've been more fortunate.
  matrix* covariance_mtx=matrix_multiplication( R_inverse,RiT);
  mimic_mtrx(covariance_mtx, covariance);
  //free everything
  matrix_free(A);
  matrix_free(I);
  matrix_free(R);
  matrix_free(R_inverse);
  matrix_free(RiT);
  matrix_free(covariance_mtx);
  vector_free(b);
  return sol;
}


int newton_with_jacobian(
	void f(const vector* x, vector* fx, matrix* J),
	vector* x,
	double epsilon,
  vector* calls,
  int which
){//this should use the newton method.
  /*where the parameters are:
    --I have to build this later -- function f: takes the input vector x, calculates analytically the vector f(x) and the Jacobian matrix of derivatives and puts them it into vector fx and matrix J;
    vector x: on input contains the starting point, on output becomes the latest approximation to the root;
    double epsilon: the accuracy goal.
  */
  int n=x->n;//number of dimensions
  matrix* Jac_matrix=matrix_alloc(n,n);
 int sweep=0;//count the # of loops executed
  double fnorm, fnorm2;
  fnorm2=epsilon+0.1; //let the while loop start
  vector* fx1=vector_alloc(n);
  vector* fx2=vector_alloc(n);
  while (fnorm2>epsilon){
  //  printf("Iteration started, vector is: %g, %g \n mtx is: \n" ,  x->data[0], x->data[1]);

    //printf("\n" );
    sweep++;
    f(x, fx1, Jac_matrix);calls->data[which]++;//creates new values
  //  print_matrix(Jac_matrix);
    vector* buffer=copy_vtr(fx1);
    matrix* R=GramSchmidt(Jac_matrix);
    vector* dx=qr_gs_solve(Jac_matrix, R, buffer);
    double lambda=1.0;
    //now scale dx with lambda
    for (int i=0;i<dx->n;i++){
      x->data[i]-=dx->data[i];//here Lambda is 1, so we only need to substract dx from x
    }
    f(x, fx2,Jac_matrix);calls->data[which]++;
    fnorm=vector_norm(fx1);
    fnorm2=vector_norm(fx2);
  /*  printf("fx1: %g %g \n",fx1->data[0], fx1->data[1] );
    printf("fx2: %g %g \n",fx2->data[0], fx2->data[1] );
    printf("fnorm2 = %g \t (1-lambda/2)*fnorm = %g \t %g is fnorm\n",fnorm2, (1-lambda/2)*fnorm , fnorm);
    */while(  fnorm2 > (1-lambda/2)*fnorm && lambda > 0.015625){
  //    printf("fnorm2 = %g \t (1-lambda/2)*fnorm = %g \t %g is fnorm\n",fnorm2, (1-lambda/2)*fnorm , fnorm);
    //  printf("%g is lambda inside\n", lambda );
      lambda=lambda/2;
      for (int i=0;i<dx->n;i++){
        dx->data[i]=dx->data[i]*1/2;
        x->data[i]-=dx->data[i];//rescaling lambda again and again...
      }
      f(x, fx2, Jac_matrix);calls->data[which]++;
      fnorm2=vector_norm(fx2);
    }
  // printf("%g is lambda\n", lambda );
    //free all the garbage created in the loop
    matrix_free(R);
    vector_free(dx);
    vector_free(buffer);
  }
  //free all the other garbage we created;
  vector_free(fx1);
  vector_free(fx2);
  matrix_free(Jac_matrix);
//  printf("fnorm 2 : %g\n", fnorm2 );
  return sweep;
}


void calc_jacobian(void f(const vector* x, vector* fx ), matrix* J, vector* x, vector* fx, double dx, vector* calls, int which){
  vector* fx_dx=vector_alloc(x->n);
  for (int i = 0; i < x->n; i++){
    x->data[i]+=dx;
    f(x, fx_dx);calls->data[which]++;
    x->data[i]-=dx;
    for (int j = 0; j < x->n; ++j) {
      double df_jdx_i = (fx_dx->data[j] - fx->data[j]) / dx;
      add_matrix_value(J, j, i, df_jdx_i);
    }
  }
  vector_free(fx_dx);
}

int newtonB(
	void f(const vector* x, vector* fx),
	vector* x,
	double epsilon,
  vector* calls,
  int which,
  double deltax
){//this should use the newton method.
  /*where the parameters are:
    --I have to build this later -- function f: takes the input vector x, calculates analytically the vector f(x) and the Jacobian matrix of derivatives and puts them it into vector fx and matrix J;
    vector x: on input contains the starting point, on output becomes the latest approximation to the root;
    double epsilon: the accuracy goal.
    SIDENOTE - I still have the same f function since I don't feel like rewriting them w/o the J matrix, I simply just neglect using J
  */
  int n=x->n;//number of dimensions
  matrix* Jac_matrix=matrix_alloc(n,n);
 int sweep=0;//count the # of loops executed
  double fnorm, fnorm2;
  fnorm2=epsilon+0.1; //let the while loop start
  vector* fx1=vector_alloc(n);
  vector* fx2=vector_alloc(n);
  while (fnorm2>epsilon){
  //  printf("Iteration started, vector is: %g, %g \n mtx is: \n" ,  x->data[0], x->data[1]);

    //printf("\n" );
    sweep++;
    f(x, fx1);calls->data[which]++;//creates new values
  //  print_matrix(Jac_matrix);
    vector* buffer=copy_vtr(fx1);
    calc_jacobian(f, Jac_matrix, x, fx1, deltax, calls, which);
    matrix* R=GramSchmidt(Jac_matrix);
    vector* dx=qr_gs_solve(Jac_matrix, R, buffer);
    double lambda=1.0;
    //now scale dx with lambda
    for (int i=0;i<dx->n;i++){
      x->data[i]-=dx->data[i];//here Lambda is 1, so we only need to substract dx from x
    }
    f(x, fx2);calls->data[which]++;
    fnorm=vector_norm(fx1);
    fnorm2=vector_norm(fx2);
    //"one should stop iterations if the step-size becomes smaller than the dx parameter"
    double smallerstep;
    if (x->data[0]>x->data[1]) smallerstep=x->data[1];
    else smallerstep=x->data[0];
  /*  printf("fx1: %g %g \n",fx1->data[0], fx1->data[1] );
    printf("fx2: %g %g \n",fx2->data[0], fx2->data[1] );
    printf("fnorm2 = %g \t (1-lambda/2)*fnorm = %g \t %g is fnorm\n",fnorm2, (1-lambda/2)*fnorm , fnorm);
    */while(  fnorm2 > (1-lambda/2)*fnorm && lambda*smallerstep > deltax){
  //    printf("fnorm2 = %g \t (1-lambda/2)*fnorm = %g \t %g is fnorm\n",fnorm2, (1-lambda/2)*fnorm , fnorm);
    //  printf("%g is lambda inside\n", lambda );
      lambda=lambda/2;
      for (int i=0;i<dx->n;i++){
        dx->data[i]=dx->data[i]*1/2;
        x->data[i]-=dx->data[i];//rescaling lambda again and again...
      }
      f(x, fx2);calls->data[which]++;
      fnorm2=vector_norm(fx2);
    }
  // printf("%g is lambda\n", lambda );
    //free all the garbage created in the loop
    matrix_free(R);
    vector_free(dx);
    vector_free(buffer);
  }
  //free all the other garbage we created;
  vector_free(fx1);
  vector_free(fx2);
  matrix_free(Jac_matrix);
//  printf("fnorm 2 : %g\n", fnorm2 );
  return sweep;
}
