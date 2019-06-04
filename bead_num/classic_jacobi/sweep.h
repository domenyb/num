#include "bases.h"
#include <math.h>


int max_of_a_row(matrix*A, int row){
  int max=0;
  double val=0;
  for (int i=row+1;i<A->cols;i++){
    double buff=extr_matrix_value(A, row, i);
    if (buff*buff>val){val=buff*buff; max=i;}
  }
  return max;
}


int Jacobi_cycles(matrix * A, vector* eigenvls, matrix* eigenvtrs ){
  int sweep_nr=0;

  int chng=1;
  for (int i=0;i<A->cols;i++){
    eigenvls->data[i]=extr_matrix_value(A,i,i);
  }

  int maxes[A->cols];
  for (int i=0;i<A->cols;i++){
    maxes[i]=max_of_a_row(A, i);
  }
  while(chng!=0){
  //build the Jpq
  chng=0;//track the changes
  sweep_nr++;
  for (int p=0;p<A->cols-1;p++)
/*    for (int q=p+1;q<A->rows;q++)*/{
    int q=maxes[p];
      double app, aqq, apq, phi;
      app=eigenvls->data[p];
      aqq=eigenvls->data[q];
      apq=extr_matrix_value(A, p, q);
     phi=0.5*atan2(/*-1**/2*apq,/*-1**/(aqq-app));//smallest
  //    else { phi=0.5*atan2(-1*2*apq,-1*(aqq-app));}//largest - the angle's orientation is changed to approach from negative direction
  //    printf("phi=%g, Apq=%g, also_ertek=%g \n ", phi, 2*apq, aqq-app);
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
    //    printf("%d %d %d %d %d\n",maxes[0],maxes[1],maxes[2],maxes[3],maxes[4] );
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
      maxes[q]=max_of_a_row(A,q);
      maxes[p]=max_of_a_row(A,p);
    }
  //}
}return sweep_nr;
  //End of function, return whether we had changes or not
}
