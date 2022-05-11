// Result_dynamicNetwork.c

// This file is used to construct series of dynamic networks from the edgelists with full temporal information

// the index of each snapshots 

#include <stdio.h>
#include <string.h>
#include "DirectedSF.h"

int process_file(int source[], int target[], int file_length) {

	char sufex[100] = "../Network_Datasets/";
	char name[100] = "science_ordered.csv";
	//"02_worrying_ordered.csv";
	//"03_zibi_ordered.csv";
	//"group16530-new-network_ordered.csv";  
	//"group151898-new-network_ordered.csv"; //"fly_vs_free_ordered.csv";
	FILE* fp = fopen(strcat(sufex,name),"r"); 
	
	if(fp==NULL){
		printf("Cannot open file\n");
		return -1;
	}

	int k;
	char buf[32];
	for(k=0; k<file_length; k++) {
		fgets(buf,sizeof(buf),fp);
		//printf("The text line %s",buf);

		int i;
		for(i = 0; i < strlen(buf); i++)
			if(buf[i] == ' ')
				break;

		int mid = i; //find the location of the ' '
		//printf("The middle character is in  %d\n",mid);

		for(i = mid+1; i < strlen(buf); i++)
			if(buf[i] == ' ')
				break;
		int mid2 = i; // for these lines with three columns
		//printf("The second middle character is in  %d\n",mid2);

		char first[16];
		char second[16];

		// handle first
		memcpy(first, buf, mid);
		first[mid] = '\0';

		int sec_size = mid2 - mid-1; 

		memcpy(second, buf + mid + 1, sec_size);
		second[sec_size] = '\0';

		//printf("The front set is [%s] \n", first);
		//printf("The behind set is [%s] \n", second);

		source[k] = atoi(first);
		//printf("The source is %d\n", source[k]);
		target[k] = atoi(second);
		//printf("The target is %d\n", target[k]);
		//printf("The num is [%d %d]\n", source[k], target[k]);

	}
	fclose(fp);
	return 0;
}

int main(void){

	char samp_name[100] = "../Network_Datasets/science_samples.txt";
	// "../Network_Datasets/02_worrying_samples.txt";
	// "../Network_Datasets/03_zibi_samples.txt";
	//"../Network_Datasets/group16530-new-network_samples.txt"; 
	//"../Network_Datasets/group151898-new-network_samples.txt"; // "../Network_Datasets/fly_vs_free_samples.txt"; 
	//char output_name[100] = "../Network_Datasets/fly_vs_free_netAttri.csv";
	FILE* indexF = fopen(samp_name,"r"); 
	//FILE* outputF = fopen(output_name,"w"); 
	
	if(indexF==NULL){
		printf("Cannot open file\n");
		return -1;
	}
	int k;
	int sample_size = 84; //92; //109; //105; //75; //108 // the number of total snapshots
	int samples[sample_size+1];

	char buf[32];
	// read in the max index of each snapshot, and build
	for (k=0; k<sample_size; k++) {
		fgets(buf,sizeof(buf),indexF);
		samples[k] = atoi(buf);
	}  

	float result[propNum];
	for (k=0; k < sample_size; k++){
		//printf("*********** Analyzing the %d-th snapshot *********** \n", k+1);
		//int* source = (int*) calloc(samples[k], sizeof(int));
		int source[samples[k]];
		int target[samples[k]];
		//int* target = (int*) calloc(samples[k], sizeof(int));
		process_file(source, target, samples[k]);
		result_t result;

		GraphAnalysis(source,target,samples[k],&result);
		//write_result(result,propNum);
		//char v[32];
		//snprintf(v,sizeof(v),"%d",result.v);
		//fputs(v, outputF);
		/*fputs(",", outputF);
		char e[32];
		snprintf(e,sizeof(e),"%d",result.e);
		fputs(e, outputF);
		fputs(",", outputF);
		char lp[32];
		snprintf(lp,sizeof(lp),"%d",result.loops);
		fputs(lp, outputF);
		fputs(",", outputF);
		fputs("\n", outputF);*/
		printf("v,e,loops,multips,maxindeg,maxoutdeg,trans,reci,reci_1n,reci_ii,reci_oo,reci_a,corr_1n,
			corr_ii,corr_io,corr_oi,corr_oo,assor,alphaIn,alphaOut,inPass,outPass,scc,in,out,tendrils,tubes,disc\n");
		printf("%d,",result.v); //The number of vertices is, 
    	printf("%d,",result.e); //The number of edges is, 
    	printf("%d,",result.loops); //The number of loops is, 
    	printf("%d,",result.multips); //The number of multiple edges is, 
    	printf("%d,", result.maxindeg); //Maximum in-degree is, 
    	printf("%d,",result.maxoutdeg); //Maximum out-degree is, 
    	printf("%f,",result.trans); //The clustering coefficient is, 
    	printf("%f,",result.reci); //The reciprocity is, 
    	printf("%f,",result.reci_1n); //The 1-node contribution is, 
    	printf("%f,",result.reci_ii); //The 2-node:i/i is is, 
    	printf("%f,",result.reci_oo); //The 2-node:o/o is is, 
    	printf("%f,",result.reci_a); //The density is, 
    	printf("%f,",result.corr_1n); //The 1-node in/out correlation is, 
    	printf("%f,",result.corr_ii); //The 2-node in/in correlation is, 
    	printf("%f,",result.corr_io); //The 2-node in/out correlation is, 
    	printf("%f,",result.corr_oi); //The 2-node out/in correlation is, 
    	printf("%f,",result.corr_oo); //The 2-node out/out correlation is,
    	printf("%f,",result.assor); //The assortativity is, 
    	printf("%.5f,",result.alphaIn); //The in-alpha is, 
    	printf("%.5f,",result.alphaOut); //The out-alpha is, 
    	printf("%d,",result.inPass); //The in-pass is, 
    	printf("%d,",result.outPass); //The out-pass is, 
    	printf("%f,",result.scc); //The scc is, 
    	printf("%f,",result.in); //The in is, 
    	printf("%f,",result.out); //The out is, 
    	printf("%f,",result.tendrils); //The tendrils is, 
    	printf("%f,",result.tubes); //The tubes is, 
   	 	printf("%f",result.disc); //The disconnect is, 
    	printf("\n");

		//printf("Done \n");
		//free(source);
		//free(target);
		//source = NULL;
		//target = NULL;
	}
	fclose(indexF);
	//fclose(outputF);
}