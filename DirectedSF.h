// DirectedSF.h

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <igraph/igraph.h>

#define propNum 28 // 20 add 4 reci + 4 corre //14 // add the bow tie results, so 14+6=20
#define repeatNum 100 // 8 // why have to be 100 and 25 ??

#define threadNum 4
// repeatPerThread = repeatNum / threadNum
#define repeatPerThread 25 // 2 //
#define initial_t 4




typedef struct result_s {
    int v;
    int e;
    int loops;
    int multips;
    int maxindeg;
    int maxoutdeg;
    int inPass;
    int outPass;
    float trans;
    // a sereis of contribution of reci
    float reci;
    float reci_1n;
    float reci_ii;
    float reci_oo;
    float reci_a;

    // a series of correlations
    float corr_1n;
    float corr_ii;
    float corr_io;
    float corr_oi;
    float corr_oo;

    float assor;
    float alphaIn;
    float alphaOut;
    // add the result of bow tie
    float scc;
    float in;
    float out;
    float tendrils;
    float tubes;
    float disc;
} result_t;

float random_num_generator();

int write_result(float result[], int length);
int write_result_lines(float result[], int length);

int count_bool_vector(igraph_vector_bool_t *v);

//float pearson(igraph_vector_t *x, igraph_vector_t *y);

//float pearson(int* a, int* b, int length); 
float pearson(int a[], int b[], int length);

int igraph_bowtie_structure(igraph_t *graph, float bowtie[]);

float count_mutual_sequence(
        int* K,
        int* Q,     
        int* seqK,
        int* seqQ,
        int v,
        int T);

int correlation_analysis(
        int* source,
        int* target,
        int* indeg,
        int* outdeg,
        int v, // number of nodes
        int T, // number of edges
        float correlation[]);

int GraphAnalysis( int* source, int* target, int T,result_t* result);

int DirectedSF(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               int* source,
               int* target);

int DirectedSF_Swap(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               int* source,
               int* target);

int DirectedSF_Local_Select(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               float phi, // phi is used to control whether select the target globally or locally
               int* source,
               int* target);

int DirectedSF_Local_Select_Swap(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               float phi, // phi is used to control whether select the target globally or locally
               int* source,
               int* target);

void* DirectedSF_thread(void* argv);
void* DirectedSF_Swap_thread(void* argv);
void* DirectedSF_Local_Select_thread(void* argv);
void* DirectedSF_Local_Select_Swap_thread(void* argv);

int DirectedSF_Result(
		int t,
		float alpha,
		float beta,
		float gama,
		float theta_in,
		float theta_out,
		float* result);

int DirectedSF_Swap_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float* result);

int DirectedSF_Local_Select_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float phi,
        float* result);

int DirectedSF_Local_Select_Swap_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float phi,
        float* result);


#if 0
unsigned long Seed=100000;
unsigned long A=48271L;           
unsigned long M=2147483647L;
unsigned int Q=44488;
unsigned int R=3399;

float random_num_generator() /** May have some problem, every time, run it, get the same random sequence **/
{
    long TmpSeed;
    TmpSeed=A*(Seed % Q)-R*(Seed /Q);
    if (TmpSeed>=0)
        Seed=TmpSeed;
    else
        Seed=TmpSeed+M;

    return (float)Seed/M;

}

void print_vector(igraph_vector_t *v, FILE *f) {
  int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

void print_bool_vector(igraph_vector_bool_t *v, FILE *f) {
  int i;
  for (i=0; i<igraph_vector_bool_size(v); i++) {
    fprintf(f, " %i", (int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}


void print_result(const igraph_plfit_result_t* result) {
  printf("====================\n");
  printf("continuous = %s\n", result->continuous ? "true" : "false");
  printf("alpha = %.5f\n", result->alpha);
  printf("xmin = %.5f\n", result->xmin);
  printf("L = %.5f\n", result->L);
  printf("D = %.5f\n", result->D);
  printf("p = %.5f\n", result->p);
  printf("====================\n");
}
#endif
