#include <stdlib.h>
int equal(double a, double b, double tau, double epsilon){
  int result=0;
  double diff=abs(a-b);
  double reldiff=diff/(abs(a)+abs(b));
  if (diff < tau) {
    result = 1;
  }
  else if (reldiff < epsilon/2) {
    result =1;
  }
  return result;
}
