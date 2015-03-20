#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "myutil.h"

void printerr(char *s)
{
	FILE *erf;
	erf = fopen("ERRFILE","w");
	printf("error:: %s\n",s);
	fprintf(erf,"error:: %s\n",s);
	exit(1);
}
