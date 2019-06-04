#include "ode.h"

double calc_func(double funs(int i, double x),int funsnr, vector* params, vector* errors, int sgn, double x){//sgn should be +1/-1/0 depending on the values
  double buffer=0;
  for (int i=0;i<funsnr;i++){
    buffer+=funs(i,x)*(params->data[i]+sgn*errors->data[i]);
  }
  return buffer;
}


int main(int argc, char** argv){
  double t=0,h=0.01, acc = 0.0, eps = 1e-3;
  vector* y = vector_alloc(2);
  FILE* f = fopen("data.dat", "w");
  for (double i = 0; i < 50; i+=0.01) {
    y->data[0]=1;
    y->data[1]=-0.5;
    t=0;//yea, I have to set them back here..............
    h=0.01;
    driver(&t, i, &h, y, acc, eps, &rkstep12, &testfunc);
          fprintf(f, "%g %g %g\n", t, y->data[0], y->data[1]);
      }
  fprintf(stdout, "The datapoints are ready for plotting. my job is done here.  \n" );
  fprintf(stdout, "A is done, let's go for B \n --------------------------------------------------------------------------------\n \n");
//innen nézd majd, most először a modified driver
  int here;
  here=0;
  t=0, acc = 0.0, eps = 1e-3;
  matrix* saved_path = matrix_alloc(10000,3);
  FILE* g = fopen("data2.dat", "w");
  //here we would start the things
  double j=50.0;
    y->data[0]=1;
    y->data[1]=-0.5;
          h = 0.01;
          t = 0;
         driver_B(&t, j, &h, y, acc, eps, &rkstep12, &testfunc,saved_path,&here);

      // At each we print the result from the matrix
      for (int i = 0; i <here ; i++) {
      fprintf(g, "%g %g %g\n", extr_matrix_value(saved_path,0, i ),  extr_matrix_value(saved_path,1,i),  extr_matrix_value(saved_path, 2,i));
      }




	vector_free(y);
	matrix_free(saved_path);
	fclose(f);
	fclose(g);
	//i don't need to free else atm.

  return 0;
}
