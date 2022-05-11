/***make the Directed_scale_free.py into Directed_scale_free.c ***/

#include "DirectedSF.h"

result_t finalResult;

// pass parameters to thread
typedef struct thread_param_s {
    int thread_id;
    int t;
    float alpha;
    float beta;
    float gama;
    float theta_in;
    float theta_out;
    float phi;
    result_t* result;
} thread_param_t;

unsigned int gSeed = 31415926;
unsigned int gStep = 3399;

// global mutex lock
pthread_mutex_t gLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t condLock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t gCond = PTHREAD_COND_INITIALIZER;


float random_num_generator(){
    unsigned ret = rand();
    srand(ret);

    return (float)rand() / (float)RAND_MAX;
}

int write_result(float result[], int length){
    // the length of the results
    if (length < propNum){
        printf("Warning, there are some results missing \n");
        return 0;
    }

    printf("The number of vertices is, %.1f \n",result[0]);
    printf("The number of edges is, %.1f \n",result[1]);
    printf("The number of loops is, %10f \n",result[2]);
    printf("The number of multiple edges is, %10f \n",result[3]);
    printf("Maximum in-degree is, %10f \n", result[4]);
    printf("Maximum out-degree is, %10f \n",result[5]);
    printf("The clustering coefficient is, %10f \n",result[6]);
    printf("The reciprocity is, %10f \n",result[7]);
    printf("The 1-node contribution is, %10f \n",result[8]);
    printf("The 2-node:i/i is is, %10f \n",result[9]);
    printf("The 2-node:o/o is is, %10f \n",result[10]);
    printf("The density is, %10f \n",result[11]);
    printf("The 1-node in/out correlation is, %10f \n",result[12]);
    printf("The 2-node in/in correlation is, %10f \n",result[13]);
    printf("The 2-node in/out correlation is, %10f \n",result[14]);
    printf("The 2-node out/in correlation is, %10f \n",result[15]);
    printf("The 2-node out/out correlation is, %10f \n",result[16]);
    printf("The assortativity is, %10f \n",result[17]);
    printf("The in-alpha is, %.5f\n",result[18]);
    printf("The out-alpha is, %.5f\n",result[19]);
    printf("The in-pass is, %.5f\n",result[20]);
    printf("The out-pass is, %.5f\n",result[21]);
    printf("The scc is, %f\n",result[22]);
    printf("The in is, %f\n",result[23]);
    printf("The out is, %f\n",result[24]);
    printf("The tendrils is, %f\n",result[25]);
    printf("The tubes is, %f\n",result[26]);
    printf("The disconnect is, %f\n",result[27]);
    printf("\n\n");

    return 0;
}

int write_result_lines(float result[], int length){
    // the length of the results
    if (length < propNum){
        printf("Warning, there are some results missing \n");
        return 0;
    }

    printf("%.1f,",result[0]);
    printf("%.1f,",result[1]);
    printf("%.1f,",result[2]);
    printf("%.1f,",result[3]);
    printf("%.1f,", result[4]);
    printf("%.1f,",result[5]);
    printf("%f,",result[6]);
    printf("%f,",result[7]);
    printf("%f,",result[8]);
    printf("%f,",result[9]);
    printf("%f,",result[10]);
    printf("%f,",result[11]);
    printf("%f,",result[12]);
    printf("%f,",result[13]);
    printf("%f,",result[14]);
    printf("%f,",result[15]);
    printf("%f,",result[16]);
    printf("%f,",result[17]);
    printf("%.5f,",result[18]);
    printf("%.5f,",result[19]);
    printf("%.1f,",result[20]);
    printf("%.1f,",result[21]);
    printf("%f,",result[22]);
    printf("%f,",result[23]);
    printf("%f,",result[24]);
    printf("%f,",result[25]);
    printf("%f,",result[26]);
    printf("%f\n",result[27]);

    return 0;
}

int count_bool_vector(igraph_vector_bool_t *v) {
  int i;
  int count=0;
  for (i=0; i<igraph_vector_bool_size(v); i++) {
    //fprintf(f, " %i", (int) VECTOR(*v)[i]);
    if(VECTOR(*v)[i]>0)
        count +=1;
  }
  //fprintf(f, "\n");
  return count;
}

int count_vector_e(int* v, int e, int length){
    // calculate if e is in v, and the number of occurences of e in v
    // if not in v, return 0; else, return the number
    int i;
    int count = 0;
    for (i=0; i < length; i++){
        //printf("Calling, the %d-th value is %d \n", i, (int) VECTOR(*v)[i] );
        if ( v[i] - e == 0)
            count += 1;
    }

    return count;
}

float pearson(int a[], int b[], int length){ // only works for sequences of integer
    // sometimes could be NAN, why?
    if (length<=0){
        printf("Warning, the length is 0\n");
    }
    double* x = (double*) calloc(length, sizeof(double));
    double* y = (double*) calloc(length, sizeof(double));

    int i;
    float pearson = 0;

    float avg_x, sum_x = 0; //average of vector x
    float avg_y, sum_y = 0;


    for(i=0; i < length; i++){
        //printf("The %d-th x is %d\n",i,a[i] );
        //x[i] = (float) a[i];
        //y[i] = (float) b[i];
        sum_x += (float) a[i];
        sum_y += (float) b[i];
    }


    avg_x = (float)sum_x/length; // igraph_vector_sum(&x)/length;// igraph_vector_size(x);
    avg_y = (float)sum_y/length;//igraph_vector_sum(&y)/length; //igraph_vector_size(y);

    for(i=0; i < length; i++){
        x[i] = (double) a[i] - avg_x;
        y[i] = (double) b[i] - avg_y;
    }

    
    double prod_0=0;
    double prod_x=0;
    double prod_y=0;

    //int i;
    for(i=0;i<length;i++){
        prod_0 += x[i] * y[i]; // igraph_vector_e(&x,i)*igraph_vector_e(&y,i);
        prod_x += x[i] * x[i];//igraph_vector_e(&x,i)*igraph_vector_e(&x,i); //(x[i]-avg_x)*(x[i]-avg_x);
        prod_y += y[i] * y[i];//igraph_vector_e(&y,i)*igraph_vector_e(&y,i);//(y[i]-avg_y)*(y[i]-avg_y);
    }
    //printf("The denominator is %f\n", prod_x*prod_y);
    if(prod_x * prod_y <= 0 ){
        //printf("Warning, the square root in pearson is 0 with length %d. \n", length);
        pearson = 0;
    }
    else if (prod_0 <= 0){
        //printf("Warning, the numerator in pearson is 0 with length %d. \n", length);

    } 
    else {
        pearson = prod_0/ sqrt(prod_x*prod_y);
    }
    

    free(x);
    free(y);


    return pearson;
}

int igraph_vector_to_array(igraph_vector_t *v, int* array, int length){
    //change from igraph_vector_t type to ordinary c array
    int i;
    for(i=0; i < length; i++){
        array[i] = (int) VECTOR(*v)[i];
    }

    return 0;
}

int igraph_bowtie_structure(igraph_t *graph, float bowtie[]){
    //return the bowtie results (6 variables) of the input graph
    //the following is to calculate clusters
    igraph_vector_t membership;
    igraph_vector_t csize;
    igraph_vector_t coreV;
    igraph_vector_t tempNeighbor;
    igraph_vector_t inV;
    igraph_vector_t outV;
    igraph_vector_t neighborIn;
    igraph_vector_t neighborOut;

    igraph_integer_t totalNodes,clusters,maxcluster,tempNode,tempSize;
    igraph_integer_t coreVnodes, inVnodes, outVnodes,neighborInnodes,neighborOutnodes,tubeVnodes=0;
    igraph_real_t scc=0, in=0, out=0, tendrils=0, tubes =0, disc=0;

    totalNodes = igraph_vcount(graph);
    //printf("The number of total nodes is %i \n", (int) totalNodes );

    igraph_vector_init(&membership,0);
    igraph_vector_init(&csize,0);
    igraph_vector_init(&coreV,0);
    igraph_vector_init(&tempNeighbor,0);
    igraph_vector_init(&inV,0);
    igraph_vector_init(&outV,0);
    igraph_vector_init(&neighborIn,0);
    igraph_vector_init(&neighborOut,0);

    igraph_clusters(graph,&membership,&csize,&clusters,IGRAPH_STRONG);
    //printf("The number of clusters is %d \n", (int)clusters);
    maxcluster = igraph_vector_which_max(&csize);
    //printf("The maximum cluster is %d \n", (int) maxcluster);
    //find the coreV nodes, whose membership equals maxcluster

    int i,j;
    for(i=0;i<igraph_vector_size(&membership);i++){
        if(VECTOR(membership)[i]==maxcluster){
            igraph_vector_push_back(&coreV,i);
        }
    }
    coreVnodes = igraph_vector_size(&coreV);
    //printf("The number of nodes in SCC is %d \n",(int) coreVnodes);
    scc = (float) coreVnodes / totalNodes;
    
    //find the neighbors of coreV, 
    for(i=0;i<coreVnodes;i++){
        //printf("The %d -th element in coreV is %d \n", i, (int)VECTOR(coreV)[i] );
        igraph_neighbors(graph,&tempNeighbor,VECTOR(coreV)[i],IGRAPH_IN);
        tempSize = igraph_vector_size(&tempNeighbor); 
        if(tempSize > 0){
            for(j=0;j<tempSize;j++){
                tempNode = VECTOR(tempNeighbor)[j];      
                if(!igraph_vector_contains(&coreV,tempNode) && !(igraph_vector_contains(&inV,tempNode))){
                    igraph_vector_push_back(&inV,tempNode);
                    //printf("The added tempNode in-neighbor is %d \n", (int) tempNode);
                }
            }
        }

        igraph_neighbors(graph,&tempNeighbor,VECTOR(coreV)[i],IGRAPH_OUT);
        tempSize = igraph_vector_size(&tempNeighbor); 
        if(tempSize > 0){
            for(j=0;j<tempSize;j++){
                tempNode = VECTOR(tempNeighbor)[j];      
                if(!igraph_vector_contains(&coreV,tempNode) && !(igraph_vector_contains(&outV,tempNode))){
                    igraph_vector_push_back(&outV,tempNode);
                    //printf("The added tempNode out-neighbor is %d \n", (int) tempNode);
                }
            }
        }

    }
    inVnodes = igraph_vector_size(&inV);
    //printf("The in-neighbor size is %d \n", (int)inVnodes);
    in = (float) inVnodes / totalNodes;

    outVnodes = igraph_vector_size(&outV);
    //printf("The out-neighbor size is %d \n", (int)outVnodes);
    out = (float) outVnodes / totalNodes; 

    if(inVnodes > 0){
        for(i=0;i< inVnodes;i++){
            //printf("The %d -th element in inV is %d \n", i, (int)VECTOR(inV)[i] );
            igraph_neighbors(graph,&tempNeighbor,VECTOR(inV)[i],IGRAPH_OUT);
            tempSize = igraph_vector_size(&tempNeighbor); 
            if(tempSize > 0){
                for(j=0;j<tempSize;j++){
                    tempNode = VECTOR(tempNeighbor)[j];      
                    if(!igraph_vector_contains(&coreV,tempNode) && !(igraph_vector_contains(&inV,tempNode))&& !(igraph_vector_contains(&outV,tempNode)))
                        igraph_vector_push_back(&neighborIn,tempNode);
                    //printf("The added tempNode in-neighbor out is %d \n", (int) tempNode);
                }
            }
        }

    }
    neighborInnodes = igraph_vector_size(&neighborIn);
    //printf("The out of in-neighbor size is %d \n", (int)neighborInnodes);


    if(outVnodes > 0){
        for(i=0;i< outVnodes;i++){
            //printf("The %d -th element in outV is %d \n", i, (int)VECTOR(outV)[i] );
            igraph_neighbors(graph,&tempNeighbor,VECTOR(outV)[i],IGRAPH_IN);
            if(igraph_vector_size(&tempNeighbor)>0){
                for(j=0;j<igraph_vector_size(&tempNeighbor);j++){
                    tempNode = VECTOR(tempNeighbor)[j];      
                    if(!igraph_vector_contains(&coreV,tempNode) && !(igraph_vector_contains(&inV,tempNode))&& !(igraph_vector_contains(&outV,tempNode)))
                        igraph_vector_push_back(&neighborOut,tempNode);
                    //printf("The added tempNode out-neighbor in is %d \n", (int) tempNode);
                }
            } 
        }
    }
    neighborOutnodes = igraph_vector_size(&neighborOut);
    //printf("The in of out-neighbor size is %d \n", (int)neighborOutnodes);
    tendrils = (float) (neighborInnodes + neighborOutnodes) / totalNodes;

    if(neighborOutnodes > 0){
        for(j=0;j< neighborOutnodes;j++){
            if(igraph_vector_contains(&neighborIn,VECTOR(neighborOut)[j])) // should be j, instead of i!!!!!!!!
                tubeVnodes += 1;
        }

    }
    //printf("The size in tubes is %d \n", tubeVnodes);
    tubes = (float) tubeVnodes / totalNodes;  

    disc = 1-scc-in-out-tendrils; // -tubes; if the tubes are contained in the tendrils
    //igraph_vector_intersect_sorted(&neighborIn,&neighborOut,&tubeV);
    //printf("The size in tubeV is %d \n", (int) igraph_vector_size(&tubeV));
    bowtie[0] = scc;
    bowtie[1] = in;
    bowtie[2] = out;
    bowtie[3] = tendrils;
    bowtie[4] = tubes;
    bowtie[5] = (disc>0) ? disc : 0;

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&csize);
    igraph_vector_destroy(&coreV);
    igraph_vector_destroy(&tempNeighbor);
    igraph_vector_destroy(&inV);
    igraph_vector_destroy(&outV);
    igraph_vector_destroy(&neighborIn);
    igraph_vector_destroy(&neighborOut);

    return 0;
}

float count_mutual_sequence(
        int* K,
        int* Q,     
        int* seqK,
        int* seqQ,
        int v,
        int T)
{
    // change to unique sequence and calculate the repeated number of each sequence
    // count the mutual links, over the number of nodes, based on the Equations
     // scan and get the unique sequence of two inputs
     int* UniK = (int*) calloc(T,sizeof(int));
     int* UniQ = (int*) calloc(T,sizeof(int));
     int* UniNum = (int*) calloc(T,sizeof(int));

     int i,j,indexUni;

     for(i=0; i < T; i++){
        UniNum[i] = 0;
     }
    
    int lengthUni = 0;
    int k_e = -1;
    int q_e = -1;

    for(i=0; i < T; i++){ 
        k_e = K[i]; 
        q_e = Q[i];
        
        indexUni = -1;
        for(j=0; j < lengthUni; j++){
            if (UniK[j] == k_e && UniQ[j] == q_e){
                indexUni = j; // the unique sequence contain the i-th sequence
            }
        }

        if (indexUni >= 0){
            UniNum[indexUni] += 1;
        }
        else {
            UniK[lengthUni] = k_e; //igraph_vector_push_back(&UniK,k_e);// append the element to the tail of the vector
            UniQ[lengthUni] = q_e; //igraph_vector_push_back(&UniQ,q_e);
            //printf("In the unique sequence, the %d-th K is %d \n", lengthUni, k_e);
            //printf("In the unique sequence, the %d-th Q is %d \n", lengthUni, q_e);
            UniNum[lengthUni] += 1;
            lengthUni += 1;
        }       
    }
    /*  
    for(i=0; i < lengthUni; i++){
        
    }
    */
    // check if i-->j is in the unique sequence, whether j-->i is also in the unique sequence 
    int invertIndex;
    float sumKQ = 0;
    
    int uk_e = -1;
    int uq_e = -1;
    int uk_node = 0;
    int uq_node = 0;

    for(i=0; i < lengthUni; i++){
        invertIndex = -1;    
        uk_e = UniK[i]; 
        uq_e = UniQ[i]; 
        //printf("The %d-th unique sequence is %d --> %d \n", i, uk_e, uq_e); 
        for(j=i+1; j < lengthUni; j++){ // exclude the self-loops
            if( (UniK[j] == uq_e) && (UniQ[j] == uk_e) ){
                invertIndex = j;
                //printf("Find the invert at %d: from %d --> %d \n",j, UniK[j], UniQ[j] );
            }
        }
        if(invertIndex >= 0){
            uk_node = count_vector_e(seqK, uk_e, v);
            uq_node = count_vector_e(seqQ, uq_e, v);
            //printf("The count of unique link is: %d \n", UniNum[i]);
            //printf("The count of invert is: %d \n", UniNum[invertIndex]);
            sumKQ = sumKQ + ((float)UniNum[i] * UniNum[invertIndex]) / (uk_node * uq_node);
            //printf("The count of node %d is: %d \n", uk_e, uk_node); 
            //printf("The count of node %d is: %d \n", uq_e, uq_node);
            //printf("The current sum is %10f \n", sumKQ);      
        }
    }

    free(UniK);
    free(UniQ);
    free(UniNum);

    return sumKQ; 
}

int correlation_analysis(
        int* source,
        int* target,       
        int* indeg,
        int* outdeg,     
        int v, // number of nodes
        int T, // number of edges
        float correlation[]) // the last one is to store the result, totally 9
{

    int* Ki = (int*) calloc(T,sizeof(int)); // the in-degree of the source node
    int* Ko = (int*) calloc(T,sizeof(int)); // the out-degree of the source node
    int* Qi = (int*) calloc(T,sizeof(int)); // the in-degree of the target node
    int* Qo = (int*) calloc(T,sizeof(int)); // the out-degree of the target node

    int i,j;
    float sumInOut = 0;
    float averInOut = 0;
    float sumKoQo = 0;
    float sumKiQi = 0;

    for(i=0; i < v; i++){
        //printf("Calling, the %d-th inseq is %d \n", i, indeg[i]);
        //printf("Calling the %d-th outseq is %d \n", i, outdeg[i]);
        sumInOut += indeg[i] * outdeg[i];
    }
    averInOut = sumInOut/v;
    //printf("The average inout is %f \n", averInOut );


    for(i=0; i < T; i++){
        Ki[i] = indeg[source[i]];
        Ko[i] = outdeg[source[i]];
        Qi[i] = indeg[target[i]];
        Qo[i] = outdeg[target[i]];
    }
    
    sumKoQo = count_mutual_sequence(Ko, Qo, outdeg, outdeg, v, T);
    sumKiQi = count_mutual_sequence(Ki, Qi, indeg, indeg, v, T);

    correlation[0] = pearson(indeg, outdeg, v);
    /*
    if (correlation[0] == 0){
        printf("The 1-node correlation is 0 \n");
        for (i = 0; i < v; ++i)
        {
            printf("The %d-th indegree is, %d\n",i,indeg[i]);          
        }
        for (i = 0; i < v; ++i){
            printf("The %d-th outdegree is, %d\n",i,outdeg[i]);
        }
    }*/

    correlation[1] = pearson(Ki, Qi, T);
    /*
    if (correlation[1] == 0){
        printf("The 2-node in:in correlation\n");
        for(i=0; i < T; i++){
            printf("The %d-th Ki is, %d \n",i, Ki[i]);
        }
        for(i=0; i < T; i++){
            printf("The %d-th Qi is, %d \n",i, Qi[i]);
        }
    }*/
    correlation[2] = pearson(Ki, Qo, T);
    /*
    if (correlation[2] == 0){
        printf("The 2-node in:out correlation\n");
        for(i=0; i < T; i++){
            printf("The %d-th Ki is, %d \n",i, Ki[i]);
        }
        for(i=0; i < T; i++){
            printf("The %d-th Qo is, %d \n",i, Qo[i]);
        }
    }*/
    correlation[3] = pearson(Ko, Qi, T);
    /*
    if (correlation[3] == 0){
        printf("The 2-node out:in correlation\n");
        for(i=0; i < T; i++){
            printf("The %d-th Ko is, %d \n",i, Ko[i]);
        }
        for(i=0; i < T; i++){
            printf("The %d-th Qi is, %d \n",i, Qi[i]);
        }
    }*/
    correlation[4] = pearson(Ko, Qo, T);
    /*
    if (correlation[4] == 0){
        printf("The 2-node out:out correlation\n");
        for(i=0; i < T; i++){
            printf("The %d-th Ko is, %d \n",i, Ko[i]);
        }
        for(i=0; i < T; i++){
            printf("The %d-th Qo is, %d \n",i, Qo[i]);
        }
    }*/
    correlation[5] = pow(averInOut,2) * pow(v,2) / pow(T,3); // r-1n
    //printf("The R1n is %lf \n", correlation[5] );
    correlation[6] = sumKoQo/T; // r 2n:o/o
    correlation[7] = sumKiQi/T; // r 2n:i/i

   // correlation[8] = (double) T /  (v*(v-1));
    //printf("The a is %f \n", correlation[8] );

    free(Ki);
    free(Ko);
    free(Qi);
    free(Qo);

    return 1;
}

int GraphAnalysis(
        int* source,// the first column in links
        int* target, // the second column in links
        int T, // the number of total links 
        result_t* result) //use to save the network properties
{
    igraph_vector_t edges;
    igraph_vector_t inseq;
    igraph_vector_t outseq;
    float bowtie[6] = {0,0,0,0,0,0};
    float corre[8] = {0,0,0,0,0,0,0,0};


    igraph_vector_bool_t loop;
    igraph_vector_bool_t multi;

    igraph_real_t trans;
    igraph_real_t density;
    igraph_real_t reci;
    igraph_real_t assor;

    igraph_plfit_result_t inpl;
    igraph_plfit_result_t outpl;

    igraph_vector_init(&edges,2*T);// length 2*T
    igraph_vector_init(&inseq,0);
    igraph_vector_init(&outseq,0);
    //igraph_vector_init(&bowtie,6);
    igraph_vector_bool_init(&loop,0);
    igraph_vector_bool_init(&multi,0);

    // initiate the vector of edges
    int i,j=0;
    for(i=0;i<T;i++){
        VECTOR(edges)[j]= (int)source[i];
        j+=1;
        VECTOR(edges)[j]= (int)target[i];
        j+=1;
    }

    igraph_t graph;
    igraph_create(&graph,&edges,0,IGRAPH_DIRECTED);
    result->v = igraph_vcount(&graph);
    result->e = igraph_ecount(&graph);

    //printf("The number of vertices is, %10i \n",(int)igraph_vcount(&graph));
    //printf("The number of edges is, %10i \n",(int)igraph_ecount(&graph));

    igraph_is_loop(&graph,&loop,igraph_ess_all(IGRAPH_EDGEORDER_ID));
    //print_bool_vector(&loop,stdout);
    //int loops = count_bool_vector(&loop);
    //printf("The number of loops is, %10f \n",loops );
    //result->loops = loops;
    result->loops = count_bool_vector(&loop);

    igraph_is_multiple(&graph,&multi,igraph_ess_all(IGRAPH_EDGEORDER_ID));
    //print_bool_vector(&multi,stdout);
    //int multips = count_bool_vector(&multi);
    //printf("The number of multiple edges is, %10i \n",multips );
    result->multips = count_bool_vector(&multi);

    //*******************************************////
    //igraph_simplify(&graph, 1, 0, 0); // This is important, it decides whether to count loops and multiple edges
    //*******************************************////

    igraph_transitivity_undirected(&graph,&trans,IGRAPH_TRANSITIVITY_NAN);
    //printf("The clustering coefficient is, %10f \n",(float)trans);
    result->trans = (float)trans;

    igraph_density(&graph,&density,0);
    result->reci_a = (float) density;

    //******!!!!!!******* This is important, whether to ignore the loops, loops may take a large portion****8!!!!
    igraph_reciprocity(&graph,&reci,1,IGRAPH_RECIPROCITY_RATIO);
    //printf("The reciprocity is, %10f \n",(double)reci);
    result->reci = (float)reci;

    igraph_degree(&graph,&inseq,igraph_vss_all(),IGRAPH_IN,1);//IGRAPH_LOOPS);
    igraph_degree(&graph,&outseq,igraph_vss_all(),IGRAPH_OUT,1);//IGRAPH_LOOPS);
    // test insequence
    //int i;

   
    int* indeg = (int*)calloc((int) result->v,sizeof(int));
    int* outdeg = (int*)calloc((int) result->v,sizeof(int));
    
    igraph_vector_to_array(&inseq, indeg, result->v);
    igraph_vector_to_array(&outseq, outdeg, result->v);


    //float pears= pearson(indeg,outdeg,result->v);
    correlation_analysis(source,target,indeg,outdeg,result->v,T,corre); // calculate the 9 correlation coefficients
    //printf("The pearson correlation between in-/out-degree is %10f \n", pears );
    //result->corr = pears;
    // a series of correlations
    result->corr_1n = corre[0];
    result->corr_ii = corre[1];
    result->corr_io = corre[2];
    result->corr_oi = corre[3];
    result->corr_oo = corre[4];

    result->reci_1n = corre[5];
    result->reci_ii = corre[6];
    result->reci_oo = corre[7];
    //result->reci_a = corre[8];
    /*
    printf("The 1-node in/out correlation is %10f \n", corre[0] );
    printf("The 2-node in/in correlation is: %10f \n", corre[1] );
    printf("The 2-node in/out correlation is: %10f \n", corre[2] );
    printf("The 2-node out/in correlation is: %10f \n", corre[3] );
    printf("The 2-node out/out correlation is: %10f \n", corre[4] );
    printf("The 1-node contribution is: %10f \n", corre[5] );
    printf("The 2-node:o/o is: %10f \n", corre[6] );
    printf("The 2-node:i/i is: %10f \n", corre[7] );
    printf("The density is: %10f \n", corre[8] );
    */

    //printf("Maximum in-degree is %10i, vertex %2i. \n", (int)igraph_vector_max(&inseq),(int)igraph_vector_which_max(&inseq));
    result->maxindeg = igraph_vector_max(&inseq);

    
    //printf("Maximum out-degree is %10i, vertex %2i. \n", (int)igraph_vector_max(&outseq),(int)igraph_vector_which_max(&outseq));
    result->maxoutdeg = igraph_vector_max(&outseq);

    igraph_assortativity_degree(&graph,&assor,1); //directed=TRUE
    //printf("The assortativity is, %10f \n",(double)assor);
    result->assor = (float)assor;

    //igraph_assortativity_degree(&graph,&assor,0);
    //printf("The 0-assortativity is, %10f \n",(double)assor);


    pthread_mutex_lock(&gLock);
    /* There is some problem here, try to decide the best alpha, start from the xmin */
    //in multi-thread, -1 is not correct, there is problem here, xmin must be greater than zero, Invalid value
    igraph_power_law_fit(&inseq,&inpl,-1,0); 
    //igraph_power_law_fit(&inseq,&inpl,0,0); //xmin must be greater than zero
    //igraph_power_law_fit(&inseq,&inpl,1,0);
    //igraph_power_law_fit(&inseq,&inpl,2,0); // alpha becomes smaller

    //igraph_power_law_fit(&inseq,&inpl,10,0); //only for t>10K, fix the xmin to 10
    //print_result(&inpl);


    igraph_power_law_fit(&outseq,&outpl,-1,0);
    //igraph_power_law_fit(&outseq,&outpl,0,0);
    //igraph_power_law_fit(&outseq,&outpl,1,0);
    //igraph_power_law_fit(&outseq,&outpl,2,0); // alpha becomes smaller
    //igraph_power_law_fit(&outseq,&outpl,10,0);// only for t>10K, fix the xmin to 10
    //print_result(&outpl);

    pthread_mutex_unlock(&gLock);

    //printf("The in-alpha is %.5f\n", inpl.alpha);
    result->alphaIn = inpl.alpha;
    result->inPass = (inpl.p<0.05) ? 0 : 1;

    //printf("The out-alpha is %.5f\n", outpl.alpha);
    result->alphaOut = outpl.alpha;
    result->outPass = (outpl.p<0.05) ? 0 : 1;

    igraph_bowtie_structure(&graph,bowtie);
    result->scc = bowtie[0];
    result->in = bowtie[1];
    result->out = bowtie[2];
    result->tendrils = bowtie[3];
    result->tubes = bowtie[4];
    result->disc = bowtie[5];
    /*
    printf("The first element in bowtie structure is %f \n", bowtie[0]);
    printf("The second element in bowtie structure is %f \n",bowtie[1]);
    printf("The third element in bowtie structure is %f \n", bowtie[2]);
    printf("The fourth element in bowtie structure is %f \n", bowtie[3]);
    printf("The fifth element in bowtie structure is %f \n", bowtie[4]);
    printf("The sixth element in bowtie structure is %f \n", bowtie[5]);
    */

    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&inseq);
    igraph_vector_destroy(&outseq);
    //igraph_vector_destroy(&bowtie);
    igraph_destroy(&graph);

    free(indeg);
    free(outdeg);
    indeg = NULL;
    outdeg = NULL;

    return 0;
}


int DirectedSF(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               int* source,
               int* target) //generate the networks, and return to soure and target
{
    //srand((unsigned)time(NULL));
    gSeed = gSeed + (unsigned)time(NULL);
    srand(gSeed);
    gSeed += gStep;

    int i = 0;
    int w = 0 ;
    int v = 0;
    int vertex = 4;
    int t0 = 4;
    int T = t + t0;
    int count = 0; //w is the target node; v is the source node

    int* indeg = (int*)calloc(T,sizeof(int));
    int* outdeg = (int*)calloc(T,sizeof(int));

    indeg[0]=indeg[1]=indeg[2]=indeg[3]=1;
    outdeg[0]=outdeg[1]=outdeg[2]=outdeg[3]=1;

    source[0]=target[3]=0;
    source[1]=target[0]=1;
    source[2]=target[1]=2;
    source[3]=target[2]=3;

    while(count < t) {
        //calculate the probability, 
        float sumIn=0;
        float sumOut=0;
        float* prop_in = (float*)calloc(vertex, sizeof(float));
        float* sump_in = (float*)calloc(vertex, sizeof(float));
        float* prop_out = (float*)calloc(vertex, sizeof(float));
        float* sump_out = (float*)calloc(vertex, sizeof(float));

        int currentT = count + t0;

        for(i = 0; i < vertex; i++) {
            prop_in[i] = (indeg[i]+theta_in)/(currentT + theta_in*vertex);
            sumIn += prop_in[i];
            sump_in[i] = sumIn;

            prop_out[i] = (outdeg[i]+theta_out)/(currentT + theta_out*vertex);
            sumOut += prop_out[i];
            sump_out[i] =sumOut;
            //printf("The in-probability and out-probability of %d vertex are %lf and %lf\n", i+1,prop_in[i],prop_out[i]);
        }

        float toss=random_num_generator();

        //printf("The random toss is %f\n", toss);
        if (toss < alpha) {
            float prob=random_num_generator();

            for(i=0;i<vertex;i++){
                if(prob <= sump_in[i])
                    break;
            }
            w = i;
            indeg[w] += 1;
            indeg[vertex] = 0;
            outdeg[vertex] = 1;

            source[currentT] = vertex;
            target[currentT] = w;
            vertex += 1;

        } else if (toss < alpha + beta) {
            float prob1=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob1 <= sump_out[i])
                    break;
            }
            v = i;

            float prob2=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob2 <= sump_in[i])
                    break;
            }
            w = i;

            outdeg[v] += 1;
            indeg[w] += 1;

            source[currentT] = v;
            target[currentT] = w;

        } else {
            float prob = random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob <= sump_out[i])
                    break;
            }
            v = i;

            outdeg[v] += 1;
            indeg[vertex] = 1;
            outdeg[vertex] = 0;

            source[currentT] = v;
            target[currentT] = vertex;

            vertex += 1;
        }

        count += 1;

        free(prop_in);
        free(sump_in);
        free(prop_out);
        free(sump_out);
    }

    free(indeg);
    free(outdeg);

    return 0;
}

int DirectedSF_Swap(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               int* source,
               int* target) //generate the networks, and return to soure and target
{
    //srand((unsigned)time(NULL));
    gSeed = gSeed + (unsigned)time(NULL);
    srand(gSeed);
    gSeed += gStep;

    int i = 0;
    int w = 0; //the target node
    int v = 0; //the source node
    int vertex = 4;
    int t0 = 4;
    int T = t + t0;
    int count = 0; //w is the target node; v is the source node

    int* indeg = (int*)calloc(T, sizeof(int));
    int* outdeg = (int*)calloc(T, sizeof(int));

    indeg[0]=indeg[1]=indeg[2]=indeg[3]=1;
    outdeg[0]=outdeg[1]=outdeg[2]=outdeg[3]=1;

    source[0]=target[3]=0;
    source[1]=target[0]=1;
    source[2]=target[1]=2;
    source[3]=target[2]=3;

    while(count<t){
        int currentT = count + t0;

        //calculate the probability, 
        float sumIn=0;
        float sumOut=0;
           float* prop_in = (float*)calloc(vertex, sizeof(float));
        float* sump_in = (float*)calloc(vertex, sizeof(float));
        float* prop_out = (float*)calloc(vertex, sizeof(float));
        float* sump_out = (float*)calloc(vertex, sizeof(float));

        for(i=0;i<vertex;i++){
            prop_in[i] = (indeg[i]+theta_in)/(currentT + theta_in*vertex);
            sumIn += prop_in[i];
            sump_in[i] = sumIn;

            prop_out[i] = (outdeg[i]+theta_out)/(currentT + theta_out*vertex);
            sumOut += prop_out[i];
            sump_out[i] =sumOut;
            //printf("The in-probability and out-probability of %d vertex are %lf and %lf\n", i+1,prop_in[i],prop_out[i]);
        }

        float toss=random_num_generator();

        if (toss < alpha){
            float prob=random_num_generator();
            //printf("In step A, the random prob  is %f\n", prob);
            for(i=0;i<vertex;i++){
                if(prob<=sump_in[i])
                    break;
            }
            w = i;
            indeg[w] += 1;
            indeg[vertex] = 0;
            outdeg[vertex] = 1;

            source[currentT] = vertex;
            target[currentT] = w;

            vertex += 1;

        } else if (toss < alpha + beta){
            float prob1=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob1<=sump_in[i]) //Original: if(prob1<=sump_out[i])
                    break;
            }
            v = i;

            float prob2=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob2<=sump_out[i]) //Original: if(prob2<=sump_in[i])
                    break;
            }
            w = i;

            outdeg[v] += 1;
            indeg[w] += 1;

            source[currentT] = v;
            target[currentT] = w;

        } else {
            float prob=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob<=sump_out[i])
                    break;
            }
            v = i;
            outdeg[v] += 1;
            indeg[vertex] = 1;
            outdeg[vertex] = 0;

            source[currentT] = v;
            target[currentT] = vertex;

            vertex += 1;
        }

        count += 1;

        free(prop_in);
        free(sump_in);
        free(prop_out);
        free(sump_out);
    }

    free(indeg);
    free(outdeg);
    indeg = NULL;
    outdeg = NULL;

    return 0;
}

int DirectedSF_Local_Select(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               float phi, // phi is used to control whether select the target globally or locally
               int* source,
               int* target) //generate the networks, and return to soure and target
{
    gSeed = gSeed + (unsigned)time(NULL);
    srand(gSeed);
    gSeed += gStep;

    int i = 0; // 
    int w = 0; // v-->w, w is the target node; 
    int v = 0; // v-->w, v is the source node
    int vertex = 4; //
    int t0 = 4; // initial number of edges
    int T = t + t0; // final number of edges
    int count = 0; // count how many edges have been added

    int* indeg = (int*)calloc(T, sizeof(int)); //save the indegree sequence in generated network
    int* outdeg = (int*)calloc(T, sizeof(int)); // 

    // initial network
    indeg[0]=indeg[1]=indeg[2]=indeg[3]=1;
    outdeg[0]=outdeg[1]=outdeg[2]=outdeg[3]=1;

    source[0]=target[3]=0;
    source[1]=target[0]=1;
    source[2]=target[1]=2;
    source[3]=target[2]=3;

    // adde one edge at one time
    while(count<t){
        int currentT = count + t0;
        //calculate the probability, 
        float sumIn=0;
        float sumOut=0;
        float* prop_in = (float*)calloc(vertex, sizeof(float));
        float* sump_in = (float*)calloc(vertex, sizeof(float));
        float* prop_out = (float*)calloc(vertex, sizeof(float));
        float* sump_out = (float*)calloc(vertex, sizeof(float));

        for(i=0;i<vertex;i++){
            prop_in[i] = (indeg[i]+theta_in)/(currentT + theta_in*vertex);
            sumIn += prop_in[i];
            sump_in[i] = sumIn;

            prop_out[i] = (outdeg[i]+theta_out)/(currentT + theta_out*vertex);
            sumOut += prop_out[i];
            sump_out[i] =sumOut;
            //printf("The in-probability and out-probability of %d vertex are %lf and %lf\n", i+1,prop_in[i],prop_out[i]);
        }

        float toss=random_num_generator();

        //printf("The random toss is %f\n", toss);
        if(toss < alpha){
            float prob=random_num_generator();
            //printf("In step A, the random prob  is %f\n", prob);
            for(i=0;i<vertex;i++){
                if(prob<=sump_in[i])
                    break;
            }
            w = i;
            indeg[w] += 1;
            indeg[vertex] = 0;
            outdeg[vertex] = 1;

            source[currentT] = vertex;
            target[currentT] = w;

            vertex += 1;

        } else if (toss < alpha + beta) {
            float prob1=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob1<=sump_out[i])
                    break;
            }
            v = i;
            outdeg[v] += 1;
            source[currentT] = v;

            //search in array, source==v or target==v,
            //remove the duplicate,
            //obtain the outdegree based on the index
            int* neighbor_v = (int*)calloc(vertex, sizeof(int)); //the maximum number of vertices
            int j=0; //index to visit the neighbors one by one until k
            int k=0; //save the current neighbor, the number of neighbors


            for(i = 0; i < currentT; i++) {
                // for each current link
                
                if(source[i] == v){
                    //printf("The source is %d\n", v );
                    for(j = 0; j < k; j++) {
                        if(neighbor_v[j] == target[i])
                            break;
                    }

                    if(j >= k){
                        neighbor_v[k] = target[i];
                        k += 1;
                        //printf("The added neighbor is target %d\n", target[i] );
                    }
                }
                else if(target[i] == v){ //only choose the one that once targets to v;
                    //printf("The target is %d\n",v);
                    for(j = 0; j < k; j++) { //check if the source has already been added into the neighbor set
                        if(neighbor_v[j] == source[i])
                            break;
                    }

                    if(j >= k){
                        neighbor_v[k] = source[i];
                        k += 1;
                        //printf("The added neighbor is source %d\n",source[i] );
                    }
                }
            }
            /*
            for(j=0;j<k;j++){
                //printf("In LSO, step B, we want to find its neighbors of  %10d \n", v);
                printf("The neighbors are: %d with indegree %d\n",neighbor_v[j],indeg[neighbor_v[j]]);
            }
            */

            //The above part is OK.

            float prob_phi = random_num_generator(); // the probability of choosing phi
                
            if((k>0) && (prob_phi<=phi)){ //at least two neighbors
                //printf("The current phi vs prob_phi is: %f vs %f \n", phi,prob_phi);    
                //printf("Select target from the neighbors####\n");
                int indegNeigbhors[k]; //the indegree of nodes in the neighbor set
                int totaldegNei=0;
                for(j=0;j<k;j++){
                    indegNeigbhors[j]=indeg[neighbor_v[j]];
                    totaldegNei += indegNeigbhors[j];
                }

                float sumIn=0;
                float prop_in_nei[k];
                float sump_in_nei[k];

                for(j=0;j<k;j++){
                    prop_in_nei[j] = (indegNeigbhors[j]+theta_in) / (totaldegNei+theta_in*k);
                    sumIn += prop_in_nei[j];
                    //printf("The %d-th neighbors' probability is %f\n", j, prop_in_nei[j]);
                    sump_in_nei[j] =sumIn;
                }

                float prob2=random_num_generator();
                for(j=0;j<k;j++){
                    if(prob2<=sump_in_nei[j])
                        break;
                }
                w = neighbor_v[j];
                //printf("The %d-th neighbor %d is chosen \n",j,w );

            }
            else{
                float prob2=random_num_generator();
                for(i=0;i<vertex;i++){
                    if(prob2<=sump_in[i])
                        break;
                }
                w = i;
            }

            free(neighbor_v);
            neighbor_v = NULL;

            indeg[w] += 1;
            //printf("In step B, the edge is from %d to %d\n", v,w);
            target[currentT] = w;
        }
        else {
            float prob=random_num_generator();
            //printf("In step C, the random prob is %f \n", prob);
            for(i=0;i<vertex;i++){
                if(prob<=sump_out[i])
                    break;
            }
            v = i;
            outdeg[v] += 1;
            indeg[vertex] = 1;
            outdeg[vertex] = 0;

            source[currentT] = v;
            target[currentT] = vertex;

            vertex += 1;
        }

        count += 1;        

        free(prop_in);
        free(sump_in);
        free(prop_out);
        free(sump_out);
    }

    free(indeg);
    free(outdeg);
    indeg = NULL;
    outdeg = NULL;

    return 0;
}

int DirectedSF_Local_Select_Swap(int t,
               float alpha,
               float beta,
               float gama,
               float theta_in,
               float theta_out,
               float phi, // phi is used to control whether select the target globally or locally
               int* source,
               int* target) //generate the networks, and return to soure and target
{
    //srand((unsigned)time(NULL));
    gSeed = gSeed + (unsigned)time(NULL);
    srand(gSeed);
    gSeed += gStep;

    int i = 0;
    int w = 0 ;
    int v = 0;
    int vertex = 4;
    int t0 = 4;
    int T = t + t0;
    int count = 0; //w is the target node; v is the source node

    int* indeg = (int*)calloc(T, sizeof(int));
    int* outdeg = (int*)calloc(T, sizeof(int));

    indeg[0]=indeg[1]=indeg[2]=indeg[3]=1;
    outdeg[0]=outdeg[1]=outdeg[2]=outdeg[3]=1;

    source[0]=target[3]=0;
    source[1]=target[0]=1;
    source[2]=target[1]=2;
    source[3]=target[2]=3;

    while(count < t){
        int currentT = count + t0;
        //calculate the probability, 
        float sumIn=0;
        float sumOut=0;
        float* prop_in = (float*)calloc(vertex, sizeof(float));
        float* sump_in = (float*)calloc(vertex, sizeof(float));
        float* prop_out = (float*)calloc(vertex, sizeof(float));
        float* sump_out = (float*)calloc(vertex, sizeof(float));

        for(i=0;i<vertex;i++){
            prop_in[i] = (indeg[i]+theta_in)/(currentT + theta_in*vertex);
            sumIn += prop_in[i];
            sump_in[i] = sumIn;

            prop_out[i] = (outdeg[i]+theta_out)/(currentT + theta_out*vertex);
            sumOut += prop_out[i];
            sump_out[i] =sumOut;
            //printf("The in-probability and out-probability of %d vertex are %lf and %lf\n", i+1,prop_in[i],prop_out[i]);
        }

        float toss=random_num_generator();

        //printf("The random toss is %f\n", toss);
        if(toss<alpha){
            float prob=random_num_generator();
            //printf("In step A, the random prob  is %f\n", prob);
            for(i=0;i<vertex;i++){
                if(prob<=sump_in[i])
                    break;
            }
            w = i;

            indeg[w] += 1;
            indeg[vertex] = 0;
            outdeg[vertex] = 1;

            source[currentT] = vertex;
            target[currentT] = w;

            vertex += 1;

        } else if(toss<alpha+beta){
            float prob1=random_num_generator();
            for(i=0;i<vertex;i++){
                if(prob1<=sump_out[i])
                    break;
            }
            v = i;
            outdeg[v] += 1;
            source[currentT] = v;

            //search in array, source==v or target==v,
            //remove the duplicate,
            //obtain the outdegree based on the index
            int* neighbor_v = (int*)calloc(vertex,sizeof(int)); //the maximum number of vertices
            int j=0;
            int k=0; //save the current neighbor
            //memset(neighbor_v, vertex+1, vertex);
            //printf("In LSO, step B, the source is %d and we want to find its neighbors \n", v );
            for(i = 0 ;i < currentT; i++){
                if(source[i] == v || target[i] == v){
                    if(source[i]==v){
                        //printf("The source is %d\n", v );
                        for(j=0;j<k;j++) {
                            if(neighbor_v[j]==target[i])
                                break;
                        }

                        if(j>=k){
                            neighbor_v[k]=target[i];
                            k += 1;
                            //printf("The added neighbor is target %d\n", target[i] );
                        }
                    }
                    else {
                        //printf("The target is %d\n",v);
                        for(j=0;j<k;j++) {
                            if(neighbor_v[j]==source[i])
                                break;
                        }

                        if(j>=k){
                            neighbor_v[k]=source[i];
                            k += 1;
                            //printf("The added neighbor is source %d\n",source[i] );
                        }
                    }
                }
            }
            /*
            for(j=0;j<k;j++){
                printf("The neighbors are: %d with outdegree %d\n",neighbor_v[j],outdeg[neighbor_v[j]]);
            }
            */            
            float prob_phi = random_num_generator(); // the probability of choosing phi
            if((k>0)&&(prob_phi<=phi)){
                //printf("Select target from neighbors####\n");
                int outdegNeigbhors[k];
                int totaldegNei=0;
                for(j=0;j<k;j++){
                    outdegNeigbhors[j]=outdeg[neighbor_v[j]];
                    totaldegNei += outdegNeigbhors[j];
                }
                float sumOut=0;
                float prop_out_nei[k];
                float sump_out_nei[k];
                for(j=0;j<k;j++){
                    prop_out_nei[j] = (outdegNeigbhors[j]+theta_out)/(totaldegNei+theta_out*k);
                    sumOut += prop_out_nei[j];
                    //printf("The %d-th neighbors' probability is %f\n", j, prop_out_nei[j]);
                    sump_out_nei[j] =sumOut;
                }

                float prob2=random_num_generator();
                for(j=0;j<k;j++){
                    if(prob2<=sump_out_nei[j])
                        break;
                }
                w = neighbor_v[j];
                //printf("The %d-th neighbor is %d \n",j,w );

            }
            else{
                float prob2=random_num_generator();
                for(i=0;i<vertex;i++){
                    if(prob2<=sump_in[i])
                        break;
                }
                w = i;
            }
            free(neighbor_v);
            neighbor_v = NULL;

            indeg[w] += 1;
            //printf("In step B, the edge is from %d to %d\n", v,w);
            target[currentT] = w;

        }
        else {
            float prob=random_num_generator();
            //printf("In step C, the random prob is %f \n", prob);
            for(i=0;i<vertex;i++){
                if(prob<=sump_out[i])
                    break;
            }
            v = i;
            outdeg[v] += 1;
            indeg[vertex] = 1;
            outdeg[vertex] = 0;
            //printf("in step C, the edge is from %d to %d\n", v,vertex);
            source[currentT] = v;
            target[currentT] = vertex;
            vertex += 1;
        }

        count += 1;

        free(prop_in);
        free(sump_in);
        free(prop_out);
        free(sump_out);
    }

    free(indeg);
    free(outdeg);
    indeg = NULL;
    outdeg = NULL;

    return 0;
}

//runs in one single thread
void* DirectedSF_thread(void* argv) {
    thread_param_t* params = (thread_param_t*) argv;
    int threadT = params->t + initial_t;

    //printf("DirectedSF_thread ID %d, T %d run...\n", params->thread_id, threadT);

    int* source = (int*) calloc(threadT, sizeof(int));
    int* target = (int*) calloc(threadT, sizeof(int));

    // loop 25 times
    int i;
    for(i = 0; i < repeatPerThread; i++) {
        DirectedSF(params->t,
                params->alpha,
                params->beta,
                params->gama,
                params->theta_in,
                params->theta_out,
                source,
                target);

        GraphAnalysis(source,
                target,
                threadT,
                params->result + i);
    }

    free(source);
    free(target);
    source = NULL;
    target = NULL;
    pthread_mutex_lock(&condLock);
    pthread_cond_signal(&gCond);
    pthread_mutex_unlock(&condLock);
    //printf("DirectedSF_thread ID: %d exit...\n", params->thread_id);

    pthread_exit(NULL);
}

void* DirectedSF_Swap_thread(void* argv) {
    thread_param_t* params = (thread_param_t*) argv;
    int threadT = params->t + initial_t;

    //printf("DirectedSF_Swap_thread ID %d, T %d run... \n",params->thread_id, threadT);
    
    int* source = (int*) calloc(threadT, sizeof(int));
    int* target = (int*) calloc(threadT, sizeof(int));

    int i;
    for(i = 0; i < repeatPerThread; i++) {
        DirectedSF_Swap(params->t,
                params->alpha,
                params->beta,
                params->gama,
                params->theta_in,
                params->theta_out,
                source,
                target);
        GraphAnalysis(source,
                target,
                threadT,
                params->result + i);
    }

    free(source);
    free(target);
    source = NULL;
    target = NULL;
    pthread_mutex_lock(&condLock);
    pthread_cond_signal(&gCond);
    pthread_mutex_unlock(&condLock);
    //printf("DirectedSF_Swap_thread ID %d exit... \n", params->thread_id);

    pthread_exit(NULL);
}

void* DirectedSF_Local_Select_thread(void* argv) {
    thread_param_t* params = (thread_param_t*) argv;
    int threadT = params->t + initial_t;

    //printf("DirectedSF_Local_Select thread ID %d, T %d run...\n", params->thread_id, threadT);

    int* source = (int*) calloc(threadT, sizeof(int));
    int* target = (int*) calloc(threadT, sizeof(int));

    // loop 25 times
    int i;
    for(i = 0; i < repeatPerThread; i++) {
        DirectedSF_Local_Select(params->t,
                params->alpha,
                params->beta,
                params->gama,
                params->theta_in,
                params->theta_out,
                params->phi,
                source,
                target);

        GraphAnalysis(source,
                target,
                threadT,
                params->result + i);
    }

    free(source);
    free(target);
    source = NULL;
    target = NULL;
    pthread_mutex_lock(&condLock);
    pthread_cond_signal(&gCond);
    pthread_mutex_unlock(&condLock);
    //printf("DirectedSF_Local_Select_thread ID: %d exit...\n", params->thread_id);

    pthread_exit(NULL);
}

void* DirectedSF_Local_Select_Swap_thread(void* argv) {
    thread_param_t* params = (thread_param_t*) argv;
    int threadT = params->t + initial_t;

    //printf("DirectedSF_Local_Select_Swap_thread ID %d, T %d run... \n",params->thread_id, threadT);
    
    int* source = (int*) calloc(threadT, sizeof(int));
    int* target = (int*) calloc(threadT, sizeof(int));

    int i;
    for(i = 0; i < repeatPerThread; i++) {
        DirectedSF_Local_Select_Swap(params->t,
                params->alpha,
                params->beta,
                params->gama,
                params->theta_in,
                params->theta_out,
                params->phi,
                source,
                target);
        GraphAnalysis(source,
                target,
                threadT,
                params->result + i);
    }

    free(source);
    free(target);
    source = NULL;
    target = NULL;
    pthread_mutex_lock(&condLock);
    pthread_cond_signal(&gCond);
    pthread_mutex_unlock(&condLock);
    //printf("DirectedSF_Local_Select_Swap_thread ID %d exit... \n", params->thread_id);

    pthread_exit(NULL);
}


int DirectedSF_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float* result) //use the generated networks to analyze and save the result.
{
    memset(result, 0, propNum);
    memset((void*)&finalResult, 0, sizeof(finalResult));
    //printf("The number of nodes is %d\n",vertex);
    //Execute the graph analysis, the input is source and target, the output is network properties
    //int t0 = 4;
    int T = t + initial_t;

    result_t resultArray[repeatNum];
    memset(resultArray, 0, sizeof(result_t) * repeatNum);

    //repeat this analysis for 100 times;
    //let each thread loops for 25 times;
    pthread_t tid[threadNum];
    thread_param_t paramsArray[threadNum] = {
        {0, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[0]},
        {1, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[25]},
        {2, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[50]},
        {3, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[75]},
    };


    // first thread
    pthread_create(&tid[0], NULL, DirectedSF_thread, (void*)(&paramsArray[0]));
    // second thread
    pthread_create(&tid[1], NULL, DirectedSF_thread, (void*)(&paramsArray[1]));
    // third thread
    pthread_create(&tid[2], NULL, DirectedSF_thread, (void*)(&paramsArray[2]));
    // forth thread
    pthread_create(&tid[3], NULL, DirectedSF_thread, (void*)(&paramsArray[3]));


    // wait all obove threads finish
    //printf("Main thread, wait all sub thread to finish...\n");
    pthread_mutex_lock(&condLock);
    pthread_cond_wait(&gCond, &condLock);
    pthread_mutex_unlock(&condLock);
 
    int thread_index;
    for (thread_index = 0; thread_index < threadNum; thread_index++) {
        pthread_join(tid[thread_index], NULL);
    }
    //printf("Main thread, all sub thread finised, got result...\n");

    int i;
    for(i=0;i<repeatNum;i++){
        finalResult.v += resultArray[i].v;
        finalResult.e += resultArray[i].e;
        finalResult.loops += resultArray[i].loops;
        finalResult.multips += resultArray[i].multips;
        finalResult.maxindeg += resultArray[i].maxindeg;
        finalResult.maxoutdeg += resultArray[i].maxoutdeg;
        finalResult.trans += resultArray[i].trans;
        // the following may only counts the networks of power-law degree distributions
        // when t is large, a lot of networks are not power-law, why??
        //
        finalResult.reci += resultArray[i].reci;
        finalResult.reci_1n += resultArray[i].reci_1n;
        finalResult.reci_ii += resultArray[i].reci_ii;
        finalResult.reci_oo += resultArray[i].reci_oo;
        finalResult.reci_a += resultArray[i].reci_a;
        finalResult.corr_1n += resultArray[i].corr_1n;
        finalResult.corr_ii += resultArray[i].corr_ii;
        finalResult.corr_io += resultArray[i].corr_io;
        finalResult.corr_oi += resultArray[i].corr_oi;
        finalResult.corr_oo += resultArray[i].corr_oo;
        finalResult.assor += resultArray[i].assor;
        
        if(resultArray[i].inPass>0) // only count the power-laws
            finalResult.alphaIn += resultArray[i].alphaIn;
        if(resultArray[i].outPass>0)
            finalResult.alphaOut += resultArray[i].alphaOut;
        //finalResult.alphaIn += resultArray[i].alphaIn;
        //finalResult.alphaOut += resultArray[i].alphaOut;
        finalResult.inPass += resultArray[i].inPass;
        finalResult.outPass += resultArray[i].outPass;
        finalResult.scc += resultArray[i].scc;
        finalResult.in += resultArray[i].in;
        finalResult.out += resultArray[i].out;
        finalResult.tendrils += resultArray[i].tendrils;
        finalResult.tubes += resultArray[i].tubes;
        finalResult.disc += resultArray[i].disc;
    }
    result[0] = (float)finalResult.v / repeatNum;
    result[1] = (float)finalResult.e / repeatNum;
    result[2] = (float)finalResult.loops / repeatNum;
    result[3] = (float)finalResult.multips / repeatNum;
    result[4] = (float)finalResult.maxindeg / repeatNum;
    result[5] = (float)finalResult.maxoutdeg / repeatNum;
    result[6] = (float)finalResult.trans / repeatNum;
    result[7] = (float)finalResult.reci / repeatNum;
    result[8] = (float)finalResult.reci_1n / repeatNum;
    result[9] = (float)finalResult.reci_ii / repeatNum;
    result[10] = (float)finalResult.reci_oo / repeatNum;
    result[11] = (float)finalResult.reci_a / repeatNum;
    result[12] = (float)finalResult.corr_1n / repeatNum;
    result[13] = (float)finalResult.corr_ii / repeatNum;
    result[14] = (float)finalResult.corr_io / repeatNum;
    result[15] = (float)finalResult.corr_oi / repeatNum;
    result[16] = (float)finalResult.corr_oo / repeatNum;
    result[17] = (float)finalResult.assor / repeatNum;
    result[18] = (float)finalResult.alphaIn / finalResult.inPass;
    result[19] = (float)finalResult.alphaOut / finalResult.outPass;
    result[20] = (float)finalResult.inPass;
    result[21] = (float)finalResult.outPass;
    result[22] = (float)finalResult.scc / repeatNum;
    result[23] = (float)finalResult.in / repeatNum;
    result[24] = (float)finalResult.out / repeatNum;
    result[25] = (float)finalResult.tendrils / repeatNum;
    result[26] = (float)finalResult.tubes / repeatNum;
    result[27] = (float)finalResult.disc / repeatNum;

    return 0;
}

int DirectedSF_Swap_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float* result) //use the generated networks to analyze and save the result.
{
    memset(result, 0, propNum);
    memset((void*)&finalResult, 0, sizeof(finalResult));
    //printf("The number of nodes is %d\n",vertex);
    int T = t + initial_t;

    result_t resultArray[repeatNum];
    memset(resultArray, 0, sizeof(result_t) * repeatNum);

    //repeat this analysis for 100 times;
    //create 4 threads, and in each thread repeat 25 times
    pthread_t tid[threadNum];
    thread_param_t paramsArray[threadNum] = {
        {0, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[0]},
        {1, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[25]},
        {2, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[50]},
        {3, t, alpha, beta, gama, theta_in, theta_out, 0, &resultArray[75]},
    };

    // create first thread
    pthread_create(&tid[0], NULL, DirectedSF_Swap_thread, (void*)(&paramsArray[0]));
    // create second thread
    pthread_create(&tid[1], NULL, DirectedSF_Swap_thread, (void*)(&paramsArray[1]));
    // create third thread
    pthread_create(&tid[2], NULL, DirectedSF_Swap_thread, (void*)(&paramsArray[2]));
    // create fourth thread
    pthread_create(&tid[3], NULL, DirectedSF_Swap_thread, (void*)(&paramsArray[3]));

    // wait for all the 4 threads to finish
    //printf("Main thread, wait all %d sub-threads to finish... \n",threadNum);
    //sleep(10);
    pthread_mutex_lock(&condLock);
    pthread_cond_wait(&gCond, &condLock);
    pthread_mutex_unlock(&condLock);

    int thread_index;
    for (thread_index = 0; thread_index < threadNum; thread_index++){
        pthread_join(tid[thread_index], NULL);
    }
    //printf("Main thread, all sub-threads finished, got result...\n");

    int i;
    for(i=0;i<repeatNum;i++){
        finalResult.v += resultArray[i].v;
        finalResult.e += resultArray[i].e;
        finalResult.loops += resultArray[i].loops;
        finalResult.multips += resultArray[i].multips;
        finalResult.maxindeg += resultArray[i].maxindeg;
        finalResult.maxoutdeg += resultArray[i].maxoutdeg;
        finalResult.trans += resultArray[i].trans;
        finalResult.reci += resultArray[i].reci;
        finalResult.reci_1n += resultArray[i].reci_1n;
        finalResult.reci_ii += resultArray[i].reci_ii;
        finalResult.reci_oo += resultArray[i].reci_oo;
        finalResult.reci_a += resultArray[i].reci_a;
        finalResult.corr_1n += resultArray[i].corr_1n;
        finalResult.corr_ii += resultArray[i].corr_ii;
        finalResult.corr_io += resultArray[i].corr_io;
        finalResult.corr_oi += resultArray[i].corr_oi;
        finalResult.corr_oo += resultArray[i].corr_oo;
        finalResult.assor += resultArray[i].assor;
        if(resultArray[i].inPass>0) // only count the power-laws
            finalResult.alphaIn += resultArray[i].alphaIn;
        if(resultArray[i].outPass>0)
            finalResult.alphaOut += resultArray[i].alphaOut;
        finalResult.inPass += resultArray[i].inPass;
        finalResult.outPass += resultArray[i].outPass;
        finalResult.scc += resultArray[i].scc;
        finalResult.in += resultArray[i].in;
        finalResult.out += resultArray[i].out;
        finalResult.tendrils += resultArray[i].tendrils;
        finalResult.tubes += resultArray[i].tubes;
        finalResult.disc += resultArray[i].disc;
    }
    result[0] = (float)finalResult.v / repeatNum;
    result[1] = (float)finalResult.e / repeatNum;
    result[2] = (float)finalResult.loops / repeatNum;
    result[3] = (float)finalResult.multips / repeatNum;
    result[4] = (float)finalResult.maxindeg / repeatNum;
    result[5] = (float)finalResult.maxoutdeg / repeatNum;
    result[6] = (float)finalResult.trans / repeatNum;
    result[7] = (float)finalResult.reci / repeatNum;
    result[8] = (float)finalResult.reci_1n / repeatNum;
    result[9] = (float)finalResult.reci_ii / repeatNum;
    result[10] = (float)finalResult.reci_oo / repeatNum;
    result[11] = (float)finalResult.reci_a / repeatNum;
    result[12] = (float)finalResult.corr_1n / repeatNum;
    result[13] = (float)finalResult.corr_ii / repeatNum;
    result[14] = (float)finalResult.corr_io / repeatNum;
    result[15] = (float)finalResult.corr_oi / repeatNum;
    result[16] = (float)finalResult.corr_oo / repeatNum;
    result[17] = (float)finalResult.assor / repeatNum;
    result[18] = (float)finalResult.alphaIn / finalResult.inPass;
    result[19] = (float)finalResult.alphaOut / finalResult.outPass;
    result[20] = (float)finalResult.inPass;
    result[21] = (float)finalResult.outPass;
    result[22] = (float)finalResult.scc / repeatNum;
    result[23] = (float)finalResult.in / repeatNum;
    result[24] = (float)finalResult.out / repeatNum;
    result[25] = (float)finalResult.tendrils / repeatNum;
    result[26] = (float)finalResult.tubes / repeatNum;
    result[27] = (float)finalResult.disc / repeatNum;
    
    return 0;
}

int DirectedSF_Local_Select_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float phi,
        float* result) //use the generated networks to analyze and save the result.
{
    memset(result, 0, propNum);
    memset((void*)&finalResult, 0, sizeof(finalResult));
    //printf("The number of nodes is %d\n",vertex);
    //Execute the graph analysis, the input is source and target, the output is network properties
    int T = t + initial_t;

    result_t resultArray[repeatNum];
    memset(resultArray, 0, sizeof(result_t) * repeatNum);

    //repeat this analysis for 100 times;
    //let each thread loops for 25 times;
    pthread_t tid[threadNum];
    thread_param_t paramsArray[threadNum] = {
        {0, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[0]},
        {1, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[25]},
        {2, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[50]},
        {3, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[75]},
    };
    // first thread
    pthread_create(&tid[0], NULL, DirectedSF_Local_Select_thread, (void*)(&paramsArray[0]));
    // second thread
    pthread_create(&tid[1], NULL, DirectedSF_Local_Select_thread, (void*)(&paramsArray[1]));
    // third thread
    pthread_create(&tid[2], NULL, DirectedSF_Local_Select_thread, (void*)(&paramsArray[2]));
    // forth thread
    pthread_create(&tid[3], NULL, DirectedSF_Local_Select_thread, (void*)(&paramsArray[3]));

    // wait all obove threads finish
    //printf("Main thread, wait all sub thread to finish...\n");
    //sleep(10);
    pthread_mutex_lock(&condLock);
    pthread_cond_wait(&gCond, &condLock);
    pthread_mutex_unlock(&condLock);
 
    int thread_index;
    for (thread_index = 0; thread_index < threadNum; thread_index++) {
        pthread_join(tid[thread_index], NULL);
    }
    //printf("Main thread, all sub thread finised, got result...\n");

    int i;
    for(i=0;i<repeatNum;i++){
        finalResult.v += resultArray[i].v;
        finalResult.e += resultArray[i].e;
        finalResult.loops += resultArray[i].loops;
        finalResult.multips += resultArray[i].multips;
        finalResult.maxindeg += resultArray[i].maxindeg;
        finalResult.maxoutdeg += resultArray[i].maxoutdeg;
        finalResult.trans += resultArray[i].trans;
        /*
        finalResult.reci += resultArray[i].reci;
        finalResult.reci_1n += resultArray[i].reci_1n;
        finalResult.reci_ii += resultArray[i].reci_ii;
        finalResult.reci_oo += resultArray[i].reci_oo;
        finalResult.reci_a += resultArray[i].reci_a;
        finalResult.corr_1n += resultArray[i].corr_1n;
        finalResult.corr_ii += resultArray[i].corr_ii;
        finalResult.corr_io += resultArray[i].corr_io;
        finalResult.corr_oi += resultArray[i].corr_oi;
        finalResult.corr_oo += resultArray[i].corr_oo;
        finalResult.assor += resultArray[i].assor;
        */
        if(resultArray[i].inPass>0) {// only count the power-laws
            finalResult.alphaIn += resultArray[i].alphaIn;
            finalResult.reci += resultArray[i].reci;
            finalResult.reci_1n += resultArray[i].reci_1n;
            finalResult.reci_ii += resultArray[i].reci_ii;
            finalResult.reci_oo += resultArray[i].reci_oo;
            finalResult.reci_a += resultArray[i].reci_a;
            finalResult.corr_1n += resultArray[i].corr_1n;
            finalResult.corr_ii += resultArray[i].corr_ii;
            finalResult.corr_io += resultArray[i].corr_io;
            finalResult.corr_oi += resultArray[i].corr_oi;
            finalResult.corr_oo += resultArray[i].corr_oo;
            finalResult.assor += resultArray[i].assor;
        }
        if(resultArray[i].outPass>0)
            finalResult.alphaOut += resultArray[i].alphaOut;
        finalResult.inPass += resultArray[i].inPass;
        finalResult.outPass += resultArray[i].outPass;
        finalResult.scc += resultArray[i].scc;
        finalResult.in += resultArray[i].in;
        finalResult.out += resultArray[i].out;
        finalResult.tendrils += resultArray[i].tendrils;
        finalResult.tubes += resultArray[i].tubes;
        finalResult.disc += resultArray[i].disc;
    }
    result[0] = (float)finalResult.v / repeatNum;
    result[1] = (float)finalResult.e / repeatNum;
    result[2] = (float)finalResult.loops / repeatNum;
    result[3] = (float)finalResult.multips / repeatNum;
    result[4] = (float)finalResult.maxindeg / repeatNum;
    result[5] = (float)finalResult.maxoutdeg / repeatNum;
    result[6] = (float)finalResult.trans / repeatNum;
    
    result[7] = (float)finalResult.reci / finalResult.inPass; //repeatNum;
    result[8] = (float)finalResult.reci_1n / finalResult.inPass;// repeatNum;
    result[9] = (float)finalResult.reci_ii / finalResult.inPass; //repeatNum;
    result[10] = (float)finalResult.reci_oo / finalResult.inPass; //repeatNum;
    result[11] = (float)finalResult.reci_a / finalResult.inPass; //repeatNum;
    result[12] = (float)finalResult.corr_1n / finalResult.inPass; //repeatNum;
    result[13] = (float)finalResult.corr_ii / finalResult.inPass; //repeatNum;
    result[14] = (float)finalResult.corr_io / finalResult.inPass; //repeatNum;
    result[15] = (float)finalResult.corr_oi / finalResult.inPass; //repeatNum;
    result[16] = (float)finalResult.corr_oo / finalResult.inPass; //repeatNum;
    result[17] = (float)finalResult.assor / finalResult.inPass; //repeatNum;

    result[18] = (float)finalResult.alphaIn / finalResult.inPass;
    result[19] = (float)finalResult.alphaOut / finalResult.outPass;
    result[20] = (float)finalResult.inPass;
    result[21] = (float)finalResult.outPass;
    result[22] = (float)finalResult.scc / repeatNum;
    result[23] = (float)finalResult.in / repeatNum;
    result[24] = (float)finalResult.out / repeatNum;
    result[25] = (float)finalResult.tendrils / repeatNum;
    result[26] = (float)finalResult.tubes / repeatNum;
    result[27] = (float)finalResult.disc / repeatNum;

    return 0;
}

int DirectedSF_Local_Select_Swap_Result(
        int t,
        float alpha,
        float beta,
        float gama,
        float theta_in,
        float theta_out,
        float phi,
        float* result) //use the generated networks to analyze and save the result.
{
    memset(result, 0, propNum);
    memset((void*)&finalResult, 0, sizeof(finalResult));
    //printf("The number of nodes is %d\n",vertex);
    //Execute the graph analysis, the input is source and target, the output is network properties
    int T = t + initial_t;

    result_t resultArray[repeatNum];
    memset(resultArray, 0, sizeof(result_t) * repeatNum);

    //repeat this analysis for 100 times;
    //let each thread loops for 25 times;
    pthread_t tid[threadNum];
    thread_param_t paramsArray[threadNum] = {
        {0, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[0]},
        {1, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[25]},
        {2, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[50]},
        {3, t, alpha, beta, gama, theta_in, theta_out, phi, &resultArray[75]},
    };
    // first thread
    pthread_create(&tid[0], NULL, DirectedSF_Local_Select_Swap_thread, (void*)(&paramsArray[0]));
    // second thread
    pthread_create(&tid[1], NULL, DirectedSF_Local_Select_Swap_thread, (void*)(&paramsArray[1]));
    // third thread
    pthread_create(&tid[2], NULL, DirectedSF_Local_Select_Swap_thread, (void*)(&paramsArray[2]));
    // forth thread
    pthread_create(&tid[3], NULL, DirectedSF_Local_Select_Swap_thread, (void*)(&paramsArray[3]));

    // wait all obove threads finish
    //printf("Main thread, wait all sub thread to finish...\n");
    //sleep(10);
    pthread_mutex_lock(&condLock);
    pthread_cond_wait(&gCond, &condLock);
    pthread_mutex_unlock(&condLock);
 
    int thread_index;
    for (thread_index = 0; thread_index < threadNum; thread_index++) {
        pthread_join(tid[thread_index], NULL);
    }
    //printf("Main thread, all sub thread finised, got result...\n");

    int i;
    for(i=0;i<repeatNum;i++){
        finalResult.v += resultArray[i].v;
        finalResult.e += resultArray[i].e;
        finalResult.loops += resultArray[i].loops;
        finalResult.multips += resultArray[i].multips;
        finalResult.maxindeg += resultArray[i].maxindeg;
        finalResult.maxoutdeg += resultArray[i].maxoutdeg;
        finalResult.trans += resultArray[i].trans;
        finalResult.reci += resultArray[i].reci;
        finalResult.reci_1n += resultArray[i].reci_1n;
        finalResult.reci_ii += resultArray[i].reci_ii;
        finalResult.reci_oo += resultArray[i].reci_oo;
        finalResult.reci_a += resultArray[i].reci_a;
        finalResult.corr_1n += resultArray[i].corr_1n;
        finalResult.corr_ii += resultArray[i].corr_ii;
        finalResult.corr_io += resultArray[i].corr_io;
        finalResult.corr_oi += resultArray[i].corr_oi;
        finalResult.corr_oo += resultArray[i].corr_oo;
        finalResult.assor += resultArray[i].assor;
        if(resultArray[i].inPass>0) // only count the power-laws
            finalResult.alphaIn += resultArray[i].alphaIn;
        if(resultArray[i].outPass>0)
            finalResult.alphaOut += resultArray[i].alphaOut;
        //finalResult.alphaIn += resultArray[i].alphaIn;
        //finalResult.alphaOut += resultArray[i].alphaOut;
        finalResult.inPass += resultArray[i].inPass;
        finalResult.outPass += resultArray[i].outPass;
        finalResult.scc += resultArray[i].scc;
        finalResult.in += resultArray[i].in;
        finalResult.out += resultArray[i].out;
        finalResult.tendrils += resultArray[i].tendrils;
        finalResult.tubes += resultArray[i].tubes;
        finalResult.disc += resultArray[i].disc;
    }
    result[0] = (float)finalResult.v / repeatNum;
    result[1] = (float)finalResult.e / repeatNum;
    result[2] = (float)finalResult.loops / repeatNum;
    result[3] = (float)finalResult.multips / repeatNum;
    result[4] = (float)finalResult.maxindeg / repeatNum;
    result[5] = (float)finalResult.maxoutdeg / repeatNum;
    result[6] = (float)finalResult.trans / repeatNum;
    result[7] = (float)finalResult.reci / repeatNum;
    result[8] = (float)finalResult.reci_1n / repeatNum;
    result[9] = (float)finalResult.reci_ii / repeatNum;
    result[10] = (float)finalResult.reci_oo / repeatNum;
    result[11] = (float)finalResult.reci_a / repeatNum;
    result[12] = (float)finalResult.corr_1n / repeatNum;
    result[13] = (float)finalResult.corr_ii / repeatNum;
    result[14] = (float)finalResult.corr_io / repeatNum;
    result[15] = (float)finalResult.corr_oi / repeatNum;
    result[16] = (float)finalResult.corr_oo / repeatNum;
    result[17] = (float)finalResult.assor / repeatNum;
    result[18] = (float)finalResult.alphaIn / finalResult.inPass;
    result[19] = (float)finalResult.alphaOut / finalResult.outPass;
    result[20] = (float)finalResult.inPass;
    result[21] = (float)finalResult.outPass;
    result[22] = (float)finalResult.scc / repeatNum;
    result[23] = (float)finalResult.in / repeatNum;
    result[24] = (float)finalResult.out / repeatNum;
    result[25] = (float)finalResult.tendrils / repeatNum;
    result[26] = (float)finalResult.tubes / repeatNum;
    result[27] = (float)finalResult.disc / repeatNum;

    return 0;
}
