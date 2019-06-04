#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <pthread.h>


struct pass_param {int a,b; double result;};

void* bar(void* param){
  struct pass_param* data = (struct pass_param*)param;
  	int a=data->a;
    int b=data->b;
    double x,y,z;
    int count=0;
    double pi;
    for(int i = 0; i<a; i++) {
      unsigned int seed = time(NULL)+i*b;
      x=rand_r(&(seed))/(double)RAND_MAX;
      y=rand_r(&(seed))/(double)RAND_MAX;
      z=x*x + y*y;
      if (z<=1) {
        count++;
      }
      }
      pi=(double)count/a*4;
  	data->result=pi;
    return 0;
}

int main(){
	int n=1e6;
	int mid=n/2;

  struct pass_param data1,data2;


  data1.a=mid;
  data2.a=mid;
  data1.b=1;
  data2.b=2;


	pthread_t traad1, traad2;
	pthread_create(&traad1,NULL,bar,(void*)&data1);
  pthread_create(&traad2,NULL,bar,(void*)&data2);
	pthread_join(traad1,NULL);
  pthread_join(traad2,NULL);

  printf("my first processor used %d datapoints: %g\n",n/2,data1.result);
  printf("my second  processor used %d datapoints: %g\n",n/2,data2.result);
  double together=(data1.result+data2.result)/2;

	printf("after %d datapoints pi is: %g\n",n,together);

return 0;
}
