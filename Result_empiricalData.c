//Result_empiricalData.c

//This file is used to calculate the network properties of the empirical data

#include <stdio.h>
#include <string.h>
#include "DirectedSF.h"

#define FILE_LENGTH 6768 //9224 //2577 //40357 // 2284546 // 876993 // 1096440 // 87627 // 140778 // 2435731 //5021410  // 7600595 //1443339 //964437 // 506550 // 948464 // 905468 
//103689 // 420045 //367662// 59835 //109396 
//40357 //   
//    the length in the *fp
 

int process_file(int* source, int* target) {

	char sufex[100] = "../Network_Datasets/";
	//printf("The char length of sufex is %lu \n", strlen(sufex) );

	// OHCs
	//char name[100] = "Links_group151898-new-network.txt";  // 40357
	
	//char name[100] = "Links_fly_vs_free.txt";   // 109396

	//char name[100] = "02_worrying_ordered.csv"; //2577
	//char name[100] = "group16530-new-network_ordered.csv"; //9224 YLD
	char name[100] = "03_zibi_ordered.csv"; // 6768 AUT

	//Web
	//char name[100] = "web-BerkStan.txt"; // 7600595

	//college message, lines of three columns
	//char name[100] = "CollegeMsg.txt"; // 59835

	//Emails
	//char name[100] = "Email-Enron.txt"; // 367662
	//char name[100] = "Email-EuAll.txt"; // 420045


	//slashdot
	//char name[100] = "Slashdot0811.txt"; //905468
	//char name[100] = "Slashdot0902.txt"; //948464

	// lines of three columns
	//char name[100] = "sx-mathoverflow.txt"; // 506550
	//char name[100] = "sx-askubuntu.txt";   // 964437
	//char name[100] = "sx-superuser.txt";    // 1443339

	// Wikipedia discussion
	//char name[100] = "WikiTalk.txt";    // 5021410
	//char name[100] = "Wiki-Vote.txt";  // 103689
	//char name[100] = "wikipedia-discussions-de/out.wikipedia-discussions-de"; // 2435731 
	//char name[100] = "wiki_talk_zh/out.wiki_talk_zh"; //2284546
	
	// Reply and discussion in threads, 
	//char name[100] = "slashdot-threads/out.slashdot-threads"; //140778
	//char name[100] = "munmun_digg_reply/out.munmun_digg_reply"; // 87627
	//char name[100] = "lkml-reply/out.lkml-reply"; // 1096440
	//char name[100] = "facebook-wosn-wall/out.facebook-wosn-wall"; //876993

	//freindship

	//printf("The char length of name is %lu \n", strlen(name) );
	printf("The file name is %s \n", name);
	FILE* fp = fopen(strcat(sufex,name),"r"); 


	if(fp==NULL){
		printf("Cannot open file\n");
		return -1;
	}

	int k;
	char buf[32];
	for(k=0; k<FILE_LENGTH; k++) {
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

		// handle second, 2 columns? 3 columns?
		//int sec_size = strlen(buf) - mid - 2;
		//int sec_size = mid2 - mid-2; // this is for lines of two columns; 
		// there are also for lines of three columns separated by ' '
		// the following can work on both two and three columns of data
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
	int T = FILE_LENGTH; //40357;
	//int source[T];
	//int target[T];
    int* source = (int*) calloc(T, sizeof(int));
    int* target = (int*) calloc(T, sizeof(int));

    process_file(source, target);
	result_t result;

	GraphAnalysis(source,target,T,&result);
	printf("The number of vertices is, %d \n",result.v);
    printf("The number of edges is, %d \n",result.e);
    printf("The number of loops is, %d \n",result.loops);
    printf("The number of multiple edges is, %d \n",result.multips);
    printf("Maximum in-degree is, %d \n", result.maxindeg);
    printf("Maximum out-degree is, %d \n",result.maxoutdeg);
    printf("The clustering coefficient is, %10f \n",result.trans);
    printf("The reciprocity is, %10f \n",result.reci);
    printf("The 1-node contribution is, %10f \n",result.reci_1n);
    printf("The 2-node:i/i is is, %10f \n",result.reci_ii);
    printf("The 2-node:o/o is is, %10f \n",result.reci_oo);
    printf("The density is, %10f \n",result.reci_a);
    printf("The 1-node in/out correlation is, %10f \n",result.corr_1n);
    printf("The 2-node in/in correlation is, %10f \n",result.corr_ii);
    printf("The 2-node in/out correlation is, %10f \n",result.corr_io);
    printf("The 2-node out/in correlation is, %10f \n",result.corr_oi);
    printf("The 2-node out/out correlation is, %10f \n",result.corr_oo);
    printf("The assortativity is, %10f \n",result.assor);
    printf("The in-alpha is, %.5f\n",result.alphaIn);
    printf("The out-alpha is, %.5f\n",result.alphaOut);
    printf("The in-pass is, %d\n",result.inPass);
    printf("The out-pass is, %d\n",result.outPass);
    printf("The scc is, %f\n",result.scc);
    printf("The in is, %f\n",result.in);
    printf("The out is, %f\n",result.out);
    printf("The tendrils is, %f\n",result.tendrils);
    printf("The tubes is, %f\n",result.tubes);
    printf("The disconnect is, %f\n",result.disc);
    printf("\n\n");

	printf("Done \n");

	free(source);
	free(target);
	source = NULL;
	target = NULL;

	return 0;
}
