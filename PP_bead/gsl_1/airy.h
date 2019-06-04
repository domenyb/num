#include<gsl/gsl_sf.h>
#include<math.h>
void airy_func(void){
	for(double x=-10;x<5;x+=0.05){
		double a=gsl_sf_airy_Ai(x,1);
		double b=gsl_sf_airy_Bi(x,1);
		fprintf(stderr,"%g %g %g\n",x,a,b);
	}
}
