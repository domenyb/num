#include<stdio.h>
#include<stdlib.h>
#include"nvector.h"

nvector* nvector_alloc(int n){
  nvector* v = malloc(sizeof(nvector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
  return v;
}

void nvector_free(nvector* v){ free(v->data); free(v);} /* v->data is identical to (*v).data */

void nvector_set(nvector* v, int i, double value){ (*v).data[i]=value; }

double nvector_get(nvector* v, int i){return (*v).data[i]; }

double nvector_dot_product(nvector* v1, nvector* v2){
  double buffer=0;
  if (v1->size!=v2->size){
    printf("problem!\n");
    return -999;
  }
  else {
    for (int i=0;i<v1->size;i++){
      buffer+=v1->data[i]*v2->data[i];
    }
    return buffer;
  }

}


int nvector_equal(nvector* v1, nvector* v2){
  if (v1->size!=v2->size){
    return 0;
  }
  else {
    for (int i=0;i<v1->size;i++){
      if(v1->data[i]!=v2->data[i]){
        return 0;
      }
    }
    return 1;
  }

}

void nvector_print(char *s, nvector * v)
{
	printf("%s", s);
	for (int i = 0; i < v->size; i++)
		printf("%9.3g ", v->data[i]);
	printf("\n");
}
