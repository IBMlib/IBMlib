/* this program checks time functions */
#include <sys/types.h>
#include <sys/times.h>
#include <stdio.h>

struct timeb {
  time_t   time;
  unsigned short millitm;
  short    timezone;
  short    dstflag;
};

extern int      ftime(struct timeb *__tp);

main()
{
  struct timeb zeit;
  int resulter, i,j, k;
  printf("hi there\n");
  for (i=0;i<20;i++)
  {
    resulter=ftime(&zeit);
    printf("%d %d %huh %hdh %hdh\n",resulter, zeit.time, zeit.millitm,
      zeit.timezone, zeit.dstflag);
    for (j=1;j<5000;j++)
      k=j*50+4/j*j;
  }
}
