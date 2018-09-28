
#include "compare.h"

#include <string>

#include <math.h>


double compare(double *a,double *b,int num,std::string legend,double tol)
{
     double t = 0.0;

     for(int i = 0; i < num ;i++)
     {
         if(fabs(a[i] - b[i]) < tol)
         {
            t += 1.0;
#ifdef COMPARE_PRINTS
            printf(" i %5d a %e b %e diff %e\n",i,a[i],b[i],fabs(a[i] - b[i]));
#endif

         }
         else
         {
#ifdef COMPARE_PRINTS
        	printf("WRONG i %5d a %e b %e diff %e\n",i,a[i],b[i],fabs(a[i] - b[i]));
#endif
         }
     }

     if(num > 0) t /= num;
#ifdef COMPARE_PRINTS
     printf("%30s %.5f\n",legend.c_str(),t);
#endif
     return t;
}

int comd(double a,double b)
{
	return (fabs(a - b) < TOLERANCE);
}
