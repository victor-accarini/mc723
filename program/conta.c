#include <stdio.h>

int main(){
	int a=1, b=2, c=4, d=10;
	
	a = a + b;
	a = a + c;
	a = a + d;
	
	if (a!=17){
	  a=2;
	} else {
	  a=1;
	}
	
	return 0;
}
