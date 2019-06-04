#include"komplex.h"
#include"stdio.h"
#define TINY 1e-6

int main(){
	komplex a = {1,2}, b = {3,4};
	printf("testing komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should   = ", R);
	komplex_print("a+b actually = ", r);

  printf("testing komplex_sub...\n");
	komplex f = komplex_sub(a,b);
	komplex F = {-2,-2};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a-b should   = ", F);
	komplex_print("a-b actually = ", f);

  printf("testing komplex_set...\n");
  komplex g={7,8};
  komplex_print("now g=",g);
  komplex_set(&g,5,6);
  komplex G = {5,6};
  printf("greal= %g",5.0);
  printf("gim= %g",6.0);
  komplex_print("g should   = ", G);
  komplex_print("g actually = ", g);


/* the following is optional */
  return 0;
}
