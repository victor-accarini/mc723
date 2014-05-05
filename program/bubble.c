#include <stdio.h>

int main(){
  int v[1000], i;
  
  for (i=0; i<1000; i++){
    v[i] = 1000-i;
  }
  
  bubble(v, 1000)
  
  return 0;
}

void bubble( int v[], int qtd )
{
  int i;
  int j;
  int aux;
  int k = qtd - 1 ;
  
  for(i = 0; i < qtd; i++)
  {
    for(j = 0; j < k; j++)
    {
      if(v[j] > v[j+1])
      {
	aux = v[j];
	v[j] = v[j+1];
	v[j+1]=aux;
      }
    }
    k--;
  }
}