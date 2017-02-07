/* problemSet1 */

#include<stdio.h>
#include<math.h>

double main()
{
     int count;
     double z, y, x;
     printf("HELLO FROM C\n");
     for(count = 1; count<11; count++)
     {
          z=count;
          x=1/z;
          y=sin(x);
          printf("%d %lf %lf\n",count,x,y);
     }
return 0;
}
