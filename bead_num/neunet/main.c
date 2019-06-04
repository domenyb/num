#include "ann.h"
#define PI 3.141592653589793238462

double gaussian(double x){
  return exp(-x*x);
}

double derivative(double x){return -2.0*x*exp(-x*x);}

double antiderivate(double x){return sqrt(PI)*erf(x)/2;}

double function_to_fit(double x){
  //let's build a function: it will be a cos*exp(x*x)*ReLu function
  if (x<0) return 0;
  else return cos(x)*exp(x);
}

double deriv_func(double x){
  if (x<0) return 0;
  else return (sin(x)*exp(x)+cos(x)*exp(x));
}

double antideriv_func(double (x)){
  if (x<0) return 0;
  else return (exp(x)*(sin(x)+ cos(x))/2);
}


int main(int argc, char** argv){//most of this function is shamelessly copied from my previous excersise and is just modified for this one

  int first_network=4;
  int second_network=8;
  ann* neunet= ann_alloc(first_network, gaussian);
  ann* neunet2=ann_alloc(second_network, gaussian);


//I probably want to put this outside as a function...
  int training=1000;
  gsl_vector* vx=gsl_vector_alloc(training);
  gsl_vector* vy=gsl_vector_alloc(training);
  double a=-3.0;
  double b=3.0;//4 is just over pi, so we'll get to see a full period
  double buffer;
  double distance_points = (b-a)/((double)training-1.00);

  for (int i=0;i<training;i++){
    buffer=a+i*distance_points;
    gsl_vector_set(vx, i, buffer);
    buffer=function_to_fit(buffer);
    gsl_vector_set(vy, i, buffer);
  }

  for (int i=0;i<first_network; i++){
    gsl_vector_set_all(neunet->data, 1);
  }

  for (int i=0;i<second_network; i++){
    gsl_vector_set_all(neunet2->data, 1);
  }

  ann_train(neunet,vx,vy);
  ann_train(neunet2,vx,vy);
  printf("we got to the end of training\n" );

  int plot=1500;
  double dist_plot=(b-a)/((double)plot-1.00);
  FILE* f = fopen("data.dat", "w");
  FILE* g= fopen("data2.dat", "w");
  FILE* h= fopen("data3.dat", "w");
  double buffer_netval, buffer_netval2, y_value;
  for (int i=0;i<plot;i++){
    buffer=a+i*dist_plot;
    buffer_netval=ann_feed_forward(neunet, buffer);
    buffer_netval2=ann_feed_forward(neunet2, buffer);
    y_value=function_to_fit(buffer);
    fprintf(f, "%g %g %g %g \n", buffer, buffer_netval, buffer_netval2, y_value );

    buffer_netval=ann_feed_forwardB(neunet, buffer, derivative);
    buffer_netval2=ann_feed_forwardB(neunet2, buffer, derivative);
    y_value=function_to_fit(buffer);
    fprintf(g, "%g %g %g %g \n", buffer, buffer_netval, buffer_netval2, y_value );

    buffer_netval=ann_feed_forwardB(neunet, buffer, antiderivate);
    buffer_netval2=ann_feed_forwardB(neunet2, buffer, antiderivate);
    y_value=function_to_fit(buffer);
    fprintf(h, "%g %g %g %g \n", buffer, buffer_netval, buffer_netval2, y_value );

  }

  printf("A and B are done, printed in separate files \n" );


  fclose(f);
  fclose(g);
  fclose(h);
  ann_free(neunet);
  ann_free(neunet2);
  gsl_vector_free(vx);
  gsl_vector_free(vy);
  return 0;
}
