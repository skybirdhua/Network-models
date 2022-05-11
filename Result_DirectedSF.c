#include <stdio.h>
#include <string.h>
#include "DirectedSF.h"

int main(int argc, char** argv) {
	if (argc < 3) {
		printf("Must specifiy the t and phi!\n");
		return -1;
	} 
	printf("Begin to run ...\n");

	int t = atoi(argv[1]);
	float phi = atof(argv[2]);
	
	// the following 3 setups is used in the growth of models
	#if 0
	//the first setup, Setup 3_1
	float alpha = 0.3;
	float beta = 0.4;
	float gama = 0.3;
	float theta_in=0.1,theta_out=0.1;
	#endif 

	#if 1
	//MDD setup, Setup 3_2
	float alpha = 0.069;
	float beta =0.862; 
	float gama = 0.069;
	float theta_in=0.1,theta_out=0.1;
	#endif

	#if 0
	//MDD setup, Setup 3_3
	float alpha = 0.069;
	float beta =0.862; 
	float gama = 0.069;
	float theta_in=4.0,theta_out=4.0;
	#endif

	// the following 3 setups are used in studying the effect of theta_in and theta_out
	#if 0
	//MDD setup theta_1
	float alpha = 0.069;
	float beta =0.862; 
	float gama = 0.069;
	float theta_in=1.5,theta_out=1.5;
	#endif

	#if 0
	//MDD setup theta_2
	float alpha = 0.069;
	float beta =0.862; 
	float gama = 0.069;
	float theta_in=3.0,theta_out=3.0;
	#endif

	#if 0
	//MDD setup theta_3
	float alpha = 0.069;
	float beta =0.862; 
	float gama = 0.069;
	float theta_in=4.5,theta_out=4.5;
	#endif
	
	float result[propNum];
	printf("*********** The current t is %d ***********  \n", t);

	printf("*** The theta_in and theta_out are  %f and %f ***  \n \n \n", theta_in, theta_out);
	
/*	
	printf("######## The Original Model: #########\n");
	DirectedSF_Result(t, alpha, beta, gama, theta_in, theta_out, result);	
	write_result(result,28);

	printf("######## The Swap-direction Model: #########\n");
	DirectedSF_Swap_Result(t, alpha, beta, gama, theta_in, theta_out, result);
	write_result(result,28);
*/
    printf("######## The LSO Model: #########\n");
	printf("The current phi is: %10f \n",phi );
	DirectedSF_Local_Select_Result(t, alpha, beta, gama, theta_in, theta_out, phi,result);
	write_result(result,28);

/*	printf("######## The LSS Model: #########\n");
	DirectedSF_Local_Select_Swap_Result(t, alpha, beta, gama, theta_in, theta_out, phi,result);	
	printf("The current phi is: %10f \n",phi );
	write_result(result,28);
*/
	printf("Done\n");
	return 0;
}
