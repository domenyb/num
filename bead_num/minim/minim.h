#include "bases.h"
/* The function is obsolete from the previous excersise
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


}*/

double func2(vector* x, vector* fx, matrix* Hesse){
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
  add_matrix_value(Hesse, 0,0,j11);
  add_matrix_value(Hesse, 0,1,j12);
  add_matrix_value(Hesse, 1,0,j21);
  add_matrix_value(Hesse, 1,1,j22);//I mean, I kinda kept most of the things, so don't really bother renaming j-s...
  return (1-x1)*(1-x1)+100*(x2-x1*x1)*(x2-x1*x1);
}

double func3(vector* x, vector* fx, matrix* jacobian){
  /*It has one local maximum at x = − 0.270845 { x=-0.270845} { x=-0.270845} and y = − 0.923039 { y=-0.923039} { y=-0.923039} where f ( x , y ) = 181.617 { f(x,y)=181.617} { f(x,y)=181.617}, and four identical local minima:

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
  add_matrix_value(jacobian, 1,1,j22);// so... jacobian = Hesse now, I just don't want to rename
  return (x1*x1+x2-11)*(x1*x1+x2-11)+(x1+x2*x2-7)*(x1+x2*x2-7);

}

/*
void func1_noJ(const vector* x, vector* fx){
  double A=10000;
  fx->data[0]=A*x->data[0]*x->data[1]-1;//  A*x*y = 1 ,
  fx->data[1]=exp(-1*x->data[0])+exp(-1*x->data[1])-1-1/A;//exp(-x) + exp(-y) = 1 + 1/A,

  //we don't need to return anything

//obsolete function
}*/

double func2_noH( vector* x, vector* fx){
  //the name of these func-s are ugly, anyways solution should be (a,a^2) where in our case a=1, so (1,1)
  double x1=x->data[0];
  double x2=x->data[1];
  //f(x,y) = (1-x)^2+100(y-x^2)^2 is the function, gradient: df/dx=-2*(1-x)+100*2*(y-x^2)*(-2x), and df/dy=200(y-x^2)
  fx->data[0]=-2*(1-x1)+100*2*(x2-x1*x1)*-2*x1;
  fx->data[1]=200*(x2-x1*x1);
  return (1-x1)*(1-x1)+100*(x2-x1*x1)*(x2-x1*x1);

}

double func3_noH(vector* x, vector* fx){
  /*It has one local maximum at x = − 0.270845 { x=-0.270845} { x=-0.270845} and y = − 0.923039 { y=-0.923039} { y=-0.923039} where f ( x , y ) = 181.617 { f(x,y)=181.617} { f(x,y)=181.617}, and four identical local minima:

  f ( 3.0 , 2.0 ) = 0.0 , { f(3.0,2.0)=0.0, } f(3.0,2.0)=0.0,
  f ( − 2.805118 , 3.131312 ) = 0.0 , { f(-2.805118,3.131312)=0.0, } f(-2.805118,3.131312)=0.0,
  f ( − 3.779310 , − 3.283186 ) = 0.0 , { f(-3.779310,-3.283186)=0.0, } f(-3.779310,-3.283186)=0.0,
  f ( 3.584428 , − 1.848126 ) = 0.0. { f(3.584428,-1.848126)=0.0. } f(3.584428,-1.848126)=0.0.
  from wikipedia. We should find exactly one of these*/
  double x1=x->data[0];
  double x2=x->data[1];
  //f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2, gradient: df/dx=2*2*(x^2+y-11)*x+2*(x+y^2-7) df/dy=2(x^2+y-11)+2*2*(x+y^2-7)*y
  fx->data[0]=4*(x1*x1 + x2 - 11)*x1 + 2*(x1 + x2*x2 - 7);
  fx->data[1]=2*(x1*x1 + x2 - 11) + 4*(x1 + x2*x2 - 7)*x2;
  return (x1*x1+x2-11)*(x1*x1+x2-11)+(x1+x2*x2-7)*(x1+x2*x2-7);

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

int newton(
double f(vector* x, vector* df, matrix* H), /* f: objective function, df: gradient, H: Hessian matrix H*/
vector* xstart, /* starting point, becomes the latest approximation to the root on exit */
double eps /* accuracy goal, on exit |gradient|<eps */
){
  int step=0;
  int dims=xstart->n;
  matrix* H=matrix_alloc(dims, dims);
  vector* gradient=vector_alloc(dims);
  double norm1=eps+0.1;//so the loop will start
  while (norm1>eps){
    double fx=f(xstart, gradient, H);
    norm1=vector_norm(gradient);
    matrix* R=GramSchmidt(H);
    vector* dx=qr_gs_solve(H, R, gradient);
    double lambda=1;
    //now scale dx with lambda
    for (int i=0;i<dx->n;i++){
      xstart->data[i]-=dx->data[i];//here Lambda is 1, so we only need to substract dx from x
    }
    double Armijo=vector_vector_product(dx, gradient);//This variable is only interesting to satisfy the armijo condition. I'm not sure whether I should recalc it or not, but I will assume the gradient doesn't change that badly
    while (f(xstart, gradient, H) > (fx+ 0.00001*lambda*Armijo) && lambda> 0.015625) {
      lambda=lambda/2;
      for (int i=0;i<dx->n;i++){
        xstart->data[i]+=dx->data[i]*lambda;
      }
    }
    f(xstart, gradient, H);//recalc
    norm1=vector_norm(gradient);//finally recalc
    matrix_free(R);
    vector_free(dx);
    step++;
  }
  matrix_free(H);
  vector_free(gradient);
  return step;
}
/*void calc_jacobian(void f(const vector* x, vector* fx ), matrix* J, vector* x, vector* fx, double dx, vector* calls, int which){
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

*/


//here comes the B
void broydens(matrix* H_inv, matrix* H_new, vector* y, vector* s, double lambda){
  int n = H_inv->rows;
  for (int i=0; i<s->n;i++){
    s->data[i]=s->data[i]*lambda*-1;
  }
  double alpha = 0.0;
  vector* object=matrix_vector_appl(H_inv, s);
  alpha= vector_vector_product(y, object);
  //printf("y: %g %g object: %g %g\n",y->data[0], y->data[1], object->data[0], object->data[1] );
  if (alpha< 0.003 && alpha> 0.003) {  //if we don't like the result, make the H_inv into eigenm.
    matrix_reset(H_inv);
    vector_free(object);
//    printf("Broyden out\n" );
  return;
  }
  vector* H_inv_y = matrix_vector_appl(H_inv, y);
  vector* pszeu = vector_alloc(n);
  for (int i=0; i<s->n;i++){
    pszeu->data[i]=s->data[i]-H_inv_y->data[i];
  }
//  printf("s  %g %g \n H_inv_y %g %g \n ",s->data[0],s->data[1], H_inv_y->data[0], H_inv_y->data[1] );

  //printf("pszeu %g %g \n object %g %g \n ",pszeu->data[0],pszeu->data[1], object->data[0], object->data[1] );
  matrix* add=vector_vector_matrix(pszeu, object);
  //print_matrix(add);
  mtx_scalar_scale(add, -1.00/alpha);
//  print_matrix(add);
  matrix_addition(H_inv, add);
  vector_free(H_inv_y);
  vector_free(pszeu);
  vector_free(object);
  matrix_free(add);
}



int qsi_newton(
double f(vector* x, vector* df), /* f: objective function, df: gradient, H: Hessian matrix H*/
vector* xstart, /* starting point, becomes the latest approximation to the root on exit */
double eps /* accuracy goal, on exit |gradient|<eps */
){
  int step=0;
  int dims=xstart->n;
  matrix* H_inv=matrix_alloc(dims, dims);//this is already an identity matrix, since this is how I allocate
  matrix* H=matrix_alloc(dims, dims);
  vector* gradient=vector_alloc(dims);
  vector* gradient_save=vector_alloc(dims);
  double norm1=eps+0.1;

  while (norm1>eps && step<100000){//just for safety measures
//    printf("x is at the start:  %g  %g\n", xstart->data[0], xstart->data[1] );
    double fx=f(xstart, gradient);
    norm1=vector_norm(gradient);
  //  printf(" this is the gradient: %g %g\n",gradient->data[0], gradient->data[1] );
/*    matrix* R=GramSchmidt(H);
    vector* dx=qr_gs_solve(H, R, gradient);*/
    vector* dx=matrix_vector_appl(H_inv, gradient);
    double lambda=1.0;
    //now scale dx with lambda
    for (int i=0;i<dx->n;i++){
      xstart->data[i]-=dx->data[i];//here Lambda is 1, so we only need to substract dx from x
    }
//    printf("armijo\n" );
    double Armijo=vector_vector_product(dx, gradient);//This variable is only interesting to satisfy the armijo condition. I'm not sure whether I should recalc it or not, but I will assume the gradient doesn't change that badly
    while (f(xstart, gradient_save) > (fx+ 0.0001*lambda*Armijo) && lambda> 0.015625/*1e-5*/) {
      lambda=lambda/2.0;
      for (int i=0;i<dx->n;i++){
        xstart->data[i]+=dx->data[i]*lambda;
      }
    }
//    printf("out the inner while\n" );
    f(xstart, gradient_save);//recalc
    norm1=vector_norm(gradient_save);//finally recalc
    //matrix_free(R);
    for (int i=0;i<dims;i++){
      gradient_save->data[i]*=-1;
      gradient->data[i]+=gradient_save->data[i];
      gradient->data[i]*=-1;
    }
//    printf("calling the Broydens, %d\n", dx->n );
    broydens(H_inv, H, gradient, dx, lambda);
//    printf("called the Broydens\n" );
      for (int i=0;i<dims; i++){
      gradient_save->data[i]=gradient->data[i];
    }
    step++;
    vector_free(dx);
    }
  matrix_free(H_inv);
  matrix_free(H);
  vector_free(gradient);
  return step;
}


double custom_func(vector*x, vector*fx){
  double A=x->data[0];
  double T=x->data[1];
  double B=x->data[2];
//I hardcode the data - i Know it is not elegant the least bit, but I kinda don't have time not to hardcode it if I can
  double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
  double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
  double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
  int N = sizeof(t)/sizeof(t[0]);
  double dA=0, dB=0, dT=0, ddA=0, ddB=0, ddT=0;//these values for the gradients
  double value=0;
  double data=0;
  double dif=0;
  for (int i=0;i<N;i++){
    value=A*exp(-t[i]/T)+B;
    data=y[i];
    dif+=(value-data)/e[i]*(value-data)/e[i];
    //the gradients analitically
    dA=exp(-t[i]/T);
    dT=A*t[i]*exp(-t[i]/T)/(T*T);
    dB=1;
    ddA+=2*(value-data)*dA*(e[i]*e[i]);
    ddB+=2*(value-data)*dB*(e[i]*e[i]);
    ddT+=2*(value-data)*dT*(e[i]*e[i]);
  }
  //update the greadients
  fx->data[0]=ddA;
  fx->data[1]=ddT;
  fx->data[2]=ddB;
  return dif;

}
