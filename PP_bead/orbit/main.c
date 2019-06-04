#include<stdio.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int log_fnc(double t,const double y[], double dydt[], void* params){
  dydt[0]=y[0]*(1-y[0]);
return GSL_SUCCESS;
}

double log_sajat(double t){
  gsl_odeiv2_system system;
  system.function=log_fnc;
  system.jacobian=NULL;
  system.dimension=1;
  system.params=NULL;
  gsl_odeiv2_driver* driver;
  double hstart=0.1;



  driver = gsl_odeiv2_driver_alloc_y_new(
    &system,
    gsl_odeiv2_step_rkf45,
    hstart,
    0.00005,
    0.00005);
    double t0=0;
    double y[]={0.5};



    gsl_odeiv2_driver_apply(driver,&t0,t,y);
    gsl_odeiv2_driver_free(driver);
    return y[0];
}
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_orbit(double t, const double y[], double dydt[], void *params){
  double eps = *(double *) params;
  dydt[0]=y[1];
  dydt[1]= 1 - y[0] + eps * y[0] * y[0];
return GSL_SUCCESS;
}

double orbit_sol(double t, double eps, double uu){
  gsl_odeiv2_system system;
  system.function=ode_orbit;
  system.jacobian=NULL;
  system.dimension=2;
  system.params= (void *) &eps;
  gsl_odeiv2_driver* driver;
  double hstart=1e-3;

  driver = gsl_odeiv2_driver_alloc_y_new(
    &system,
    gsl_odeiv2_step_rkf45,
    hstart,
    0.000005,
    0.000005);
    double t0=0;
    double y[2]={1,uu};


    gsl_odeiv2_driver_apply(driver,&t0,t,y);
    gsl_odeiv2_driver_free(driver);
    return y[0];
}




int main(){
		for (double x = (0); x < 3; x += 0.1)
			fprintf(stdout, "%g %g %g\n", x, log_sajat(x), exp(x)/(exp(x)+1));
		double eps_1 = 0;
		double uu_1 = 0;
		double eps_2 = 0;
		double uu_2 = -0.5;
		double eps_3 = 0.01;
		double uu_3 = -0.5;
		for (double x = (0); x < 20 * M_PI; x += 0.05)
			fprintf(stderr, "%g %g %g %g\n", x, orbit_sol(x,eps_1,uu_1),orbit_sol(x,eps_2,uu_2),orbit_sol(x,eps_3,uu_3));


return 0;
}
