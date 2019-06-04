#include "bases.h"
#include <math.h>

int jacobi_sweep(matrix * A, vector* eigenvls, matrix* eigenvtrs){
  //build the Jpq
  int chng=0;//track the changes
  for (int p=0;p<A->rows;p++){
    for (int q=p+1;q<A->rows;q++){
      double app, aqq, apq;
      app=eigenvls->data[p];
      aqq=eigenvls->data[q];
      apq=extr_matrix_value(A, p, q);
      double phi=0.5*atan2(2*apq,aqq-app);
      double sine= sin(phi);
      double cosine= cos(phi);//I calculate these preamptively, since otherwise I'd need to call these multiple times
      double app1=cosine*cosine*app-2*sine*cosine*apq+sine*sine*aqq;
      double aqq1=sine*sine*app+2*sine*cosine*apq+cosine*cosine*aqq;
      //if don't need changes, give it back as it is
      if (app1!=app || aqq!=aqq1){
        chng=1;
        eigenvls->data[p]=app1;
        eigenvls->data[q]=aqq1;
        add_matrix_value(A, p, q, 0);
        //now let's change the values
        for (int i=0;i<p;i++){
          double aip=extr_matrix_value(A, i, p);
          double aiq=extr_matrix_value(A, i, q);
          add_matrix_value(A, i, p, cosine*aip-sine*aiq );
          add_matrix_value(A,i,q,cosine*aiq+sine*aip);
        }
        for(int i=p+1;i<q;i++){
          double api=extr_matrix_value(A,p,i);
          double aiq=extr_matrix_value(A,i,q);
          add_matrix_value(A,p,i,cosine*api-sine*aiq);
          add_matrix_value(A,i,q,cosine*aiq+sine*api);
        }
        for(int i=q+1;i<A->cols;i++){
          double api=extr_matrix_value(A,p,i);
          double aqi=extr_matrix_value(A,q,i);
          add_matrix_value(A,p,i,cosine*api-sine*aqi);
          add_matrix_value(A,q,i,cosine*aqi+sine*api);
        }
        //also in the eigenvtrs
        for(int i=0;i<A->cols;i++){
          double evtrip=extr_matrix_value(eigenvtrs,i,p);
          double evtriq=extr_matrix_value(eigenvtrs,i,q);
          add_matrix_value(eigenvtrs,i,p,cosine*evtrip-sine*evtriq);
          add_matrix_value(eigenvtrs,i,q,cosine*evtriq+sine*evtrip);
        }
      }
    }
  }
  return chng;
  //End of function, return whether we had changes or not
}
