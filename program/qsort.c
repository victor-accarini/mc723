#include <stdio.h>
#include <stdlib.h>

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

int main()
{
   int n, i;
   int values[1000];
   scanf("%d", &n);
   
   for( i = 0 ; i < n; i++ ) {
      scanf("%d ", values[i]);
   }
   
   printf("Before sorting the list is: \n");
   for( i = 0 ; i < n; i++ ) {
      printf("%d ", values[i]);
   }

   qsort(values, n, sizeof(int), cmpfunc);

   printf("\nAfter sorting the list is: \n");
   for( i = 0 ; i < n; i++ ) {
      printf("%d ", values[i]);
   }
  
  return(0);
}