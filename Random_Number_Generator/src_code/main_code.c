#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

	//N is the number of numbers to be generate
	//hist_bins is the number of bins to build the histogram
	int N=300000;
	srand(99);
	int i,total;
	float PI=3.14159265359;
	float r[N],inv[N],arg[N],bin[N];
	int hist_bins = 150,hist_b[N];

	FILE *outo;
	outo=fopen("outo.dat","w");
	fprintf(outo,"i,unif,inv,hist_b\n");

	//loop over bins
	for(i=0;i<hist_bins;i++){
		hist_b[i] = 0;	
	}	

	//applying inverse transform method and allocating numbers at bins
	for(i=0;i<N;i++){
		r[i]=(float)rand()/RAND_MAX;
		arg[i]=1-(2*r[i]);
		inv[i]=acosf(arg[i]);
		hist_b[(int)floor(inv[i]*hist_bins/PI)] +=1;
		
	}
	
	//printing info at file
	for(i=0;i<hist_bins;i++){
		fprintf(outo,"%d,%f,%f,%d\n",i, r[i], inv[i], hist_b[i]);
	}
}
