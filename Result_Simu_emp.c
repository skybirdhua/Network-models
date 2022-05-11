#include <stdio.h>
#include <string.h>
#include "DirectedSF.h"

int main(int argc, char** argv) {
	if (argc < 2) {
		printf("Must specifiy the phi!\n");
		return -1;
	} 
	printf("Begin to run ...\n");

	float phi = atof(argv[1]);

	// the following two setups are used in simulating the real-world networks

	#if 0
	float theta_in=0.1;
	float theta_out=0.1;
	#endif

	#if 1
	float theta_in=4.0;
	float theta_out=4.0;
	#endif

	// the setups in empirical data
	#if 0   //MDD setup with loops,   40357
	int t = 40353;
	float alpha = 0.0625;
	float beta =0.875; 
	float gama = 0.0625;	
	#endif

	#if 0   //MDD setup with 27-th snapshot,   272
	int t = 268;
	float alpha = 0.1838;
	float beta = 0.6324; 
	float gama = 0.1838;	
	#endif

	#if 0   //MDD setup with 40-th snapshot,   2789
	int t = 2785;
	float alpha = 0.1153;
	float beta =0.7694;
	float gama = 0.1153;	
	#endif

	#if 0   //MDD setup with 60-th snapshot,   22571
	int t = 22567;
	float alpha = 0.0682; 
	float beta =0.8636;
	float gama = 0.0682;	
	#endif

	#if 0   // Depression setup with loops,  109396
	int t = 109392; 
	float alpha = 0.0533;
	float beta = 0.8934;
	float gama = 0.0533;
	#endif

	#if 0   // Depression setup with 20-th snapshot,  289
	int t = 285; 
	float alpha = 0.1488;
	float beta = 0.7024; 
	float gama = 0.1488;
	#endif

	#if 0   // Depression setup with 57-th snapshot,  12804
	int t = 12800; 
	float alpha = 0.1053; 
	float beta = 0.7894;
	float gama = 0.1053;
	#endif	

	#if 0   // Depression setup with 70-th snapshot,  25829
	int t = 25825; 
	float alpha = 0.0775;
	float beta = 0.8450;
	float gama = 0.0775;
	#endif

	#if 0   // Depression setup with 90-th snapshot,  77494
	int t = 77490; 
	float alpha = 0.0563;
	float beta = 0.8874;
	float gama = 0.0563;
	#endif

	#if 0   //YLD setup with loops,   9224
	int t = 9220;
	float alpha = 0.1771;
	float beta =0.6457; 
	float gama = 0.1771;	
	#endif

	#if 1   //AUT setup with loops,   6768
	int t = 6764;
	float alpha = 0.1713;
	float beta =0.6574; 
	float gama = 0.1713;	
	#endif




	#if 0   // Email-EuAll,  420045
	int t = 420041;
	float alpha = 0.316;
	float beta = 0.368;
	float gama = 0.316;
	#endif

	#if 0   // CollegeMsg, edges  59835
	int t = 59831; 
	float alpha = 0.0159;
	float beta = 0.9682;
	float gama = 0.0159;
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
	
	printf("######## The LSS Model: #########\n");
	DirectedSF_Local_Select_Swap_Result(t, alpha, beta, gama, theta_in, theta_out, phi,result);	
	printf("The current phi is: %10f \n",phi );
	write_result(result,28);
	
	printf("Done\n");
	return 0;
}
