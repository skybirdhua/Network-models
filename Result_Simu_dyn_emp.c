#include <stdio.h>
#include <string.h>
#include "DirectedSF.h"

// simulate a series of dynamic networks, with an array of v and e

// write the results separately, one Model in one file

//

int main(void) {

	// the following two setups are used in simulating the real-world networks

	#if 0
	float theta_in=0.1;
	float theta_out=0.1;
	#endif

	#if 1
	float theta_in=4.0;
	float theta_out=4.0;
	#endif

	#if 0
	int e[55] = {10182,10713,11677,12804,13840,14532,15080,15999,16534,17239,17769,18288,19203,19894,20813,21697,25829,28940,32955,35638,37176,38621,40025,41687,43170,44937,48158,51164,55296,57793,60931,64165,66757,69746,72153,75038,77494,80348,82104,83793,85638,87682,89493,91652,94026,96333,98923,101514,103410,104700,105920,106868,107777,109030,109396};
	int v[55] = {2265,2370,2572,2697,2796,2889,2961,3075,3135,3207,3278,3340,3468,3547,3658,3754,4004,4198,4376,4509,4683,4832,4986,5161,5354,5533,5775,6106,6452,6778,7124,7495,7770,8057,8323,8533,8724,8931,9091,9278,9467,9702,9893,10088,10295,10522,10700,10863,11009,11128,11221,11364,11494,11624,11660};
	float phi = 0.8;
	#endif

	#if 1
	
	int e[27] ={10689,11754,12666,13871,14822,16309,17398,18605,19446,20344,21540,22571,23636,24542,25537,26733,27885,29066,30283,31358,32310,33941,35998,37437,38549,39737,40357};
	int v[27] ={1612,1721,1827,1989,2130,2298,2432,2574,2681,2810,2951,3080,3200,3318,3437,3559,3698,3828,3996,4129,4250,4431,4602,4743,4852,4973,5051};
	float phi = 0.8;
	#endif
	// the setups in empirical data

	float alpha = 0; 
	float beta = 0; 
	float gama = 0;


	printf("v,e,loops,multips,maxindeg,maxoutdeg,trans,reci,reci_1n,reci_ii,reci_oo,reci_a,corr_1n,corr_ii,corr_io,corr_oi,corr_oo,assor,alphaIn,alphaOut,inPass,outPass,scc,in,out,tendrils,tubes,disc\n");

	int i=0, t=0, length = 27;

	float result[propNum];

	
	for (i=0; i<length; i++){
		alpha = (float) v[i]/(2.0*e[i]);
		beta = 1.0 - 2.0*alpha;
		gama = alpha;
		t = e[i] - 4;

		//DirectedSF_Result(t, alpha, beta, gama, theta_in, theta_out,result);
		//DirectedSF_Swap_Result(t, alpha, beta, gama, theta_in, theta_out,result);
		//DirectedSF_Local_Select_Result(t, alpha, beta, gama, theta_in, theta_out,phi,result);
		DirectedSF_Local_Select_Swap_Result(t, alpha, beta, gama, theta_in, theta_out, phi,result);
		write_result_lines(result,28);
	}
	 

	#if 0
	for (i=0; i<length; i++){
		
		write_result_lines(result,28)
	}
	#endif 

	#if 0
	for (i=0; i<length; i++){
		
		write_result_lines(result,28)
	}
	#endif 

	#if 0
	for (i=0; i<length; i++){
		write_result_lines(result,28)
	}
	#endif 

	return 0;
}
