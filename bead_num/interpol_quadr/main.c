#include <stdio.h>
#include <stdlib.h>

typedef struct {int n/*,last*/; double *x, *y, *b, *c; } qspline; //last: the last binary search's return value, if the next Z is in the same intervall, you don't need to run the binary search again. At the moment it is not implemented

qspline * qspline_alloc(int n, double *x, double *y){
  qspline * object=malloc(sizeof(qspline));
  object->n=n;
//  object->last=0;
  object->x=malloc(sizeof(double)*n);
  object->y=malloc(sizeof(double)*n);
  object->b=malloc(sizeof(double)*(n-1));
  object->c=malloc(sizeof(double)*(n-1));
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
  object->c[0]=0;
  for (int i=0;i<n-2;i++){
    object->c[i+1]=mered[i+1]-mered[i]-object->c[i]*(period[i]/period[i+1]);
  }
  object->c[n-2]=object->c[n-2]/2;
  for (int i=n-3;i>=0;i--){
    object->c[i]=mered[i+1]-mered[i]-object->c[i]*(period[i+1]/period[i]);
  }
  for (int i=0; i<n-1;i++){
    object->b[i]=mered[i]-object->c[i]*period[i];
  }
  return object;
};
//here we created a qspline object

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

double qspline_interp(qspline * q, double z){
  double result;
  int lower_ind=binary(q->x, q->n, z);
  double subperiod=z-q->x[lower_ind];
  result=q->y[lower_ind]+subperiod*(q->b[lower_ind]+subperiod*q->c[lower_ind]);
  return result;
};
//interpolates with the given qspline

double qspline_derivative(qspline * q, double z){
  //qspline's derivative at point Z
  double result;
  int lower_ind=binary(q->x, q->n, z);
  result=q->b[lower_ind]+2*q->c[lower_ind]*z-2*q->c[lower_ind]*q->x[lower_ind];
  return result;
}

double qspline_integr(qspline * q, double z){
  double result=0;
  int lower_ind=binary(q->x, q->n, z);
  for (int i=0;i<lower_ind;i++){
    //integrate for the intervals below
    result+=(q->y[i]-q->b[i]*q->x[i]+q->c[i]*q->x[i]*q->x[i])*(q->x[i+1]-q->x[i]); //the constants integrated
    result+=(q->b[i]-2*q->c[i]*q->x[i])*1/2*(q->x[i+1]*q->x[i+1]-q->x[i]*q->x[i]);//the first order integrated
    result+=(q->c[i])*1/3*(q->x[i+1]*q->x[i+1]*q->x[i+1]-q->x[i]*q->x[i]*q->x[i]);//second order integrated
  }
  //now we've integrated all the intervals below, comes the interval "left in half"
  result+=(q->y[lower_ind]-q->b[lower_ind]*q->x[lower_ind]+q->c[lower_ind]*q->x[lower_ind]*q->x[lower_ind])*(z-q->x[lower_ind]); //the constants integrated
  result+=(q->b[lower_ind]-2*q->c[lower_ind]*q->x[lower_ind])*1/2*(z*z-q->x[lower_ind]*q->x[lower_ind]);//the first order integrated
  result+=(q->c[lower_ind])*1/3*(z*z*z-q->x[lower_ind]*q->x[lower_ind]*q->x[lower_ind]);//second order integrated
  return result;
}

void qspline_free (qspline * q){
  free(q->x);
  free(q->y);
  free(q->b);
  free(q->c);
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
  qspline * my_spline=qspline_alloc(n, x, y);
  free(x);
  free(y);
  for (int i=0;i<n;i++){
    fprintf(stderr,"%lg \t",x[i]);
    fprintf(stderr,"%lg \t",y[i]);
    fprintf(stderr,"%lg \n", qspline_integr(my_spline, x[i]));
  }
  double z,zy,zint;
  for (int i=0;i<nz;i++){
      fscanf(f, "%lg", &z);
  //    printf("%g \n \n ", z);
      zy=qspline_interp(my_spline, z);
      zint=qspline_integr(my_spline,z);
      fprintf(stdout, "%lg \t %lg \t %lg \n", z, zy, zint);

  }
  qspline_free(my_spline);

  return 0;
}
