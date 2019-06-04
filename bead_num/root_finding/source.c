unction_linear_equation_with_j(gsl_vector* x,gsl_vector* y, gsl_matrix* J, gsl_vector* functioncall){ // The system of linear equation A*x*y = 1 , exp(-x) + exp(-y) = 1 + 1/A, with A = 10000 is implemented
gsl_vector_set(functioncall,0,gsl_vector_get(functioncall,0)+1);
    double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    double A = 10000;
		gsl_vector_set(y,0, A * x_1*x_2 - 1); // We Calculate the value of the function A*x*y -1 = 0 in the given point
		gsl_vector_set(y,1, exp(-x_1) + exp(-x_2) - 1 - 1.0/A); // We Calculate the value of the function exp(-x) + exp(-y) -( 1 + 1/A)= 0  in the given point

    double J_11 = A * x_2 ; // Analytically Calculated d_X f1(x,y)
    double J_12 = A * x_1; // Analytically Calculated d_y  f1(x,y)
    double J_21 = -exp(-x_1); // Analytically Calculated d_X  f2(x,y)
    double J_22 = -exp(-x_2); // Analytically Calculated d_y f2(x,y)
    gsl_matrix_set(J, 0, 0, J_11);
    gsl_matrix_set(J, 0, 1, J_12);
    gsl_matrix_set(J, 1, 0, J_21);
    gsl_matrix_set(J, 1, 1, J_22);
}


void function_himmel_with_j(gsl_vector* x,gsl_vector* y, gsl_matrix* J, gsl_vector* functioncall){ // The Himmelblau function f(x,y) = (x2+y-11)2+(x+y2-7)2 is implemented
gsl_vector_set(functioncall,2,gsl_vector_get(functioncall,2)+1);
		double x_1=gsl_vector_get(x,0), x_2=gsl_vector_get(x,1);
    gsl_vector_set(y,0,  4*(x_1*x_1 + x_2 - 11)*x_1 + 2*(x_1 + x_2*x_2 - 7)); // We Calculate the value of the function in the given point
		gsl_vector_set(y,1, 2*(x_1*x_1 + x_2 - 11) + 4*(x_1 + x_2*x_2 - 7)*x_2); // We calculate the gradient of the function

    double J_11 = 12*x_1*x_1 + 4*x_2-42; // Analytically Calculated d_X d_X f(x,y)
    double J_12 = 4*(x_1+x_2); // Analytically Calculated d_X d_y f(x,y)
    double J_21 = 4*(x_1+x_2); // Analytically Calculated d_y d_X f(x,y)
    double J_22 = 4*x_1+12*x_2*x_2-26; // Analytically Calculated d_y d_y f(x,y)
    gsl_matrix_set(J, 0, 0, J_11);
    gsl_matrix_set(J, 0, 1, J_12);
    gsl_matrix_set(J, 1, 0, J_21);
    gsl_matrix_set(J, 1, 1, J_22);
}






nt newton_with_jacobian(void function(gsl_vector* x, gsl_vector* fx, gsl_matrix* J, gsl_vector* functioncall), gsl_vector* x, double epsilon,  gsl_vector* functioncall){
    //  This function calculates the root by Newton's method with a analytically jacobian
    // Takes in f: takes the input vector x, calculates analytically the vector f(x), the jacobian matrix, secong arguement is vector x: on input contains the starting point, on output becomes the latest approximation to the root;
    //  at last it takes double epsilon: the accuracy goal.

// The function is based on the python script example in the lecture notes:
// We allocate the needed parameters
int number_of_eq=x->size;
gsl_matrix* R = gsl_matrix_alloc(number_of_eq, number_of_eq);
gsl_matrix* J = gsl_matrix_alloc(number_of_eq,number_of_eq);
gsl_vector* fx = gsl_vector_alloc(number_of_eq);

gsl_vector* deltax = gsl_vector_alloc(number_of_eq);
gsl_vector* new_fx = gsl_vector_alloc(number_of_eq);
gsl_vector* cfx = gsl_vector_alloc(number_of_eq);
double lambda = 1.00;
int step=0;
double functionnorm, newfunctionnorm;
// We start the stepping procedure:
do {
  step++;
function(x,fx,J,functioncall); // We update the Jacobian and funciton values.
gsl_vector_memcpy(cfx, fx);
// Now we need to Solve system J∆x=−f(x) for stepping size delta x this is done through QR decomposition and multiplying with Q-invese as done in previous exercise on linear equations.
// First we QR decompose
 qr_gs_decomp(J, R);
 // Then we solve the system using backsubstitution
 qr_gs_solve(J, R, cfx, deltax);
 lambda = 1.00;
 // We find the new x by x + lambda* -deltax
gsl_vector_scale (deltax, -lambda);
gsl_vector_add(x, deltax);

// we evaluate the function at the new x
function(x,new_fx,J,functioncall);
// We find the norm of the function, being the distance to 0.0 for both the old and new function
functionnorm = gsl_blas_dnrm2(fx);
newfunctionnorm = gsl_blas_dnrm2(new_fx);

// We find the propper lambda factor by evluating the criteria
 while (newfunctionnorm > (1-lambda/2.00)*functionnorm && lambda > 1.00/64.00) {
   // If the lambda is to large, we half it and calculate new values
   lambda /= 2.00;
   gsl_vector_scale(deltax, lambda);
   gsl_vector_add(x, deltax);
   // we evaluate the function at the new x
   function(x,new_fx,J,functioncall);
   newfunctionnorm = gsl_blas_dnrm2(new_fx);
 }

// When lambda is rescaled we find the new value of the function at the new place
function(x,new_fx,J,functioncall);
newfunctionnorm = gsl_blas_dnrm2(new_fx);
}while(newfunctionnorm > epsilon);

// We free the parameters
gsl_matrix_free(R);
gsl_matrix_free(J);
gsl_vector_free(fx);
gsl_vector_free(deltax);
gsl_vector_free(cfx);
gsl_vector_free(new_fx);
return step;
}
