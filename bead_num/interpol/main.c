#include <stdio.h>
#include <stdlib.h>

int binary(double *tomb, int len, double z){
  //tomb is the vector of x or y values, len is the number of values, z is the point to which I want to interpolate. I'll assume that x and y are ordered lists.
  int lower=0, upper=len-1;
  while (upper-lower>1){
    int middle=(upper+lower)/2;
    if (tomb[middle]>z)
      upper=middle;
    else
      lower=middle;
  };
  return lower;
};


double linterp(int n, double *x, double *y, double z){
  double result;
  int lower_ind=binary(x, n, z);
  result=y[lower_ind]+(y[lower_ind+1]-y[lower_ind])/(x[lower_ind+1]-x[lower_ind])*(z-x[lower_ind]);
  return result;

};

double linterp_integ(int n, double *x, double *y, double z){
  double result=0;
  int lower_ind=binary(x,n,z);
  for (int i=0;i<lower_ind;i++){
    result+=(y[i+1]+y[i])/2*(x[i+1]-x[i]); //I realise this is not a literal integration, but for a linear spline, it is the same
//    fprintf(stderr,"result for %lg: %lg \n", z, result);
  };
  if (z!=x[lower_ind]){
    result+=(linterp(n, x, y, z)+y[lower_ind])/2*(z-x[lower_ind]);
  }
  return result;
};




int main (){
  FILE* f;//our datapoints in a file
  f= fopen("datapoints.dat", "r");
  int n, nz;//n --> nr of datapoints, nz --> number of z-s
  //  fprintf(stderr, "%i \n", n);
  fscanf(f, "%d", &n); //the file's first row contains the number of datapoints and the number of "z"-s to interpolate to
  fscanf(f, "%d", &nz);
  double* y=malloc(sizeof(double)*n);
  double* x=malloc(sizeof(double)*n);
  for (int i=0;i<n;i++){
    fscanf(f,"%lg",&x[i]);
    fscanf(f,"%lg",&y[i]);
  }
  for (int i=0;i<n;i++){
    fprintf(stderr,"%lg \t",x[i]);
    fprintf(stderr,"%lg \t",y[i]);
    fprintf(stderr,"%lg \n",linterp_integ(n, x, y, x[i]));
  }
  double z,zy,zint;
  for (int i=0;i<nz;i++){
      fscanf(f, "%lg", &z);
  //    printf("%g \n \n ", z);
      zy=linterp(n, x, y, z);
      zint=linterp_integ(n,x,y,z);
      fprintf(stdout, "%lg \t %lg \t %lg \n", z, zy, zint);

  }
  free(x);
  free(y);
  fclose(f);
  return 0;
}
