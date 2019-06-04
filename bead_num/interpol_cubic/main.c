#include <stdio.h>
#include <stdlib.h>

typedef struct {int n/*,last*/; double *x, *y, *b, *c,*d; } cspline; //last: the last binary search's return value, if the next Z is in the same intervall, you don't need to run the binary search again. At the moment it is not implemented

cspline * cspline_alloc(int n, double *x, double *y){
  cspline * object=malloc(sizeof(cspline));
  object->n=n;
//  object->last=0;
  object->x=malloc(sizeof(double)*n);
  object->y=malloc(sizeof(double)*n);
  object->b=malloc(sizeof(double)*(n-1));
  object->c=malloc(sizeof(double)*(n-1));
  object->d=malloc(sizeof(double)*(n-1));
  double period[n-1];
  double mered[n-1];
  for (int i=0;i<n-1;i++){
    period[i]=x[i+1]-x[i];
    mered[i]=(y[i+1]-y[i])/period[i];
  }
  for (int i=0;i<n;i++){
    object->x[i]=x[i];
    object->y[i]=y[i];
  };
//build the tridiagonal sys:
  double D[n], Q[n-1], B[n];
  D[0]=2;
  D[n-1]=2;
  Q[0]=1;
  B[0]=3*mered[0];
  B[n-1]=3*mered[n-2];
  for (int i=0;i<n-2;i++){
    D[i+1]=2*(period[i]/period[i+1]+1);
    Q[i+1]=period[i]/period[i+1];
    B[i+1]=3*(mered[i]+mered[i+1]*period[i]/period[i+1]);
  }
  //values estabilished, now comes the gauss elimination
  for (int i=1;i<n;i++){
    D[i]-=Q[i-1]/D[i-1];
    B[i]-=B[i-1]/D[i-1];
  }
  //and back-substitution
  object->b[n-1]=B[n-1]/D[n-1];
  for (int i=n-2;i>=0;i--) {
    object->b[i]=(B[i]-Q[i]*object->b[i+1])/D[i];
  }
  for (int i=0;i<n-1;i++){
  object->c[i]=(-2*object->b[i]-object->b[i+1]+3*mered[i])/period[i];
  object->d[i]=(object->b[i]+object->b[i+1]-2*mered[i])/period[i]/period[i];
  }
  return object;
};
//here we created a cspline object

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
//still a binary search

double cspline_interp(cspline * q, double z){
  double result;
  int lower_ind=binary(q->x, q->n, z);
  double subperiod=z-q->x[lower_ind];
  result=q->y[lower_ind]+subperiod*(q->b[lower_ind]+subperiod*(q->c[lower_ind]+subperiod*q->d[lower_ind]));
  return result;
};
//interpolates with the given cspline

double cspline_derivative(cspline * q, double z){
  //cspline's derivative at point Z
  double result;
  int lower_ind=binary(q->x, q->n, z);
  result=q->b[lower_ind]+2*q->c[lower_ind]*z-2*q->c[lower_ind]*q->x[lower_ind]+3*q->d[lower_ind]*(z-q->x[lower_ind])*(z-q->x[lower_ind]);
  return result;
}

double cspline_integr(cspline * q, double z){
  double result=0;
  int lower_ind=binary(q->x, q->n, z);
  for (int i=0;i<lower_ind;i++){
    //integrate for the intervals below
    result+=(q->y[i]-q->b[i]*q->x[i]+q->c[i]*q->x[i]*q->x[i]-q->d[i]*q->x[i]*q->x[i]*q->x[i])*(q->x[i+1]-q->x[i]); //the constants integrated
    result+=(q->b[i]-2*q->c[i]*q->x[i]+3*q->d[i]*q->x[i]*q->x[i])*1/2*(q->x[i+1]*q->x[i+1]-q->x[i]*q->x[i]);//the first order integrated
    result+=(q->c[i]-3*q->d[i]*q->x[i])*1/3*(q->x[i+1]*q->x[i+1]*q->x[i+1]-q->x[i]*q->x[i]*q->x[i]);//second order integrated
    result+=(q->d[i])*1/4*(q->x[i+1]*q->x[i+1]*q->x[i+1]*q->x[i+1]-q->x[i]*q->x[i]*q->x[i]*q->x[i]);//third order integrated
  }
  //now we've integrated all the intervals below, comes the interval "left in half"
  result+=(q->y[lower_ind]-q->b[lower_ind]*q->x[lower_ind]+q->c[lower_ind]*q->x[lower_ind]*q->x[lower_ind]-q->d[lower_ind]*q->x[lower_ind]*q->x[lower_ind]*q->x[lower_ind])*(z-q->x[lower_ind]); //the constants integrated
  result+=(q->b[lower_ind]-2*q->c[lower_ind]*q->x[lower_ind]+3*q->d[lower_ind]*q->x[lower_ind]*q->x[lower_ind])*1/2*(z*z-q->x[lower_ind]*q->x[lower_ind]);//the first order integrated
  result+=(q->c[lower_ind]-3*q->d[lower_ind]*q->x[lower_ind])*1/3*(z*z*z-q->x[lower_ind]*q->x[lower_ind]*q->x[lower_ind]);//second order integrated
  result+=(q->d[lower_ind])*1/4*(z*z*z*z-q->x[lower_ind]*q->x[lower_ind]*q->x[lower_ind]*q->x[lower_ind]);//third order integrated
  return result;
}

void cspline_free (cspline * q){
  free(q->x);
  free(q->y);
  free(q->b);
  free(q->c);
  free(q->d);
  free(q);
}


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
  cspline * my_spline=cspline_alloc(n, x, y);
  free(x);
  free(y);
  for (int i=0;i<n;i++){
    fprintf(stderr,"%lg \t",x[i]);
    fprintf(stderr,"%lg \t",y[i]);
    fprintf(stderr,"%lg \n", cspline_integr(my_spline, x[i]));
  }
  double z,zy,zint;
  for (int i=0;i<nz;i++){
      fscanf(f, "%lg", &z);
  //    printf("%g \n \n ", z);
      zy=cspline_interp(my_spline, z);
      zint=cspline_integr(my_spline,z);
      fprintf(stdout, "%lg \t %lg \t %lg \n", z, zy, zint);

  }
  cspline_free(my_spline);
  return 0;
}
