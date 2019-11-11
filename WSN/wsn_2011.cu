#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>


/***all macros**/
#define ALPHA1 50                   // nano joule 1e-9
#define ALPHA2 50                   // nano joule 1e-9
#define BETA   100                  // pico joule = 1e-12
#define EXP_M  4
#define E_INIT 5                    // in Joules

#define  MAX_NODE  100000000
#define  DEBUG 1




#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/**all type declaration***/

using namespace std;

class Node{
public:
    unsigned int val;
    vector<unsigned int> weights;
    vector<Node*> Edges;

    Node(int val){
        this->val = val;
    }

    void addEdge(Node* v,unsigned int w){
        this->Edges.push_back(v);
        this->weights.push_back(w);
    }

};

// __device__ volatile int Cx[MAX_NODE];
__device__ volatile int PQ[MAX_NODE];


//K in parallel
__global__ void extractMin(int* PQ_size, int* expandNodes,int* expandNodes_size,float* Cx,float* ln,int* openList,int N,int K){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id<K && PQ_size[id]>0){

        //extract min from PQ
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];

        // restructure the heap
        PQ[front]=PQ[front+PQ_size[id]-1];
        PQ_size[id]-=1;
        int pqIndex = 0;

        while(2*pqIndex+1 < PQ_size[id]){
            if(2*pqIndex+2 >= PQ_size[id]){            

                if( ln[PQ[front+pqIndex]] > ln[PQ[front+2*pqIndex+1]]){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else if(ln[PQ[front+pqIndex]] == ln[PQ[front+2*pqIndex+1]]){
                    
                    if( Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+1]]){
                        int swap = PQ[front + 2*pqIndex+1];
                        PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                        PQ[front + pqIndex] = swap;
                        pqIndex = 2*pqIndex+1;
                    }
                    else
                        break;
                }
                else 
                    break;
    
            }
            else{
                if(ln[PQ[front+pqIndex]] > ln[PQ[front+2*pqIndex+1]] && ln[PQ[front+2*pqIndex+1]] <= ln[PQ[front+2*pqIndex+2]] ){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else if(ln[PQ[front+pqIndex]] > ln[PQ[front+2*pqIndex+2]] && ln[PQ[front+2*pqIndex+2]] <= ln[PQ[front+2*pqIndex+1]] ){
                    int swap = PQ[front + 2*pqIndex+2];
                    PQ[front + 2*pqIndex+2] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+2;                    
                } 
                else if(ln[PQ[front+pqIndex]] == ln[PQ[front+2*pqIndex+1]] && ln[PQ[front+pqIndex]] == ln[PQ[front+2*pqIndex+2]]  ){
                   
                    if( Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+1]] && Cx[PQ[front+2*pqIndex+1]] <= Cx[PQ[front+2*pqIndex+2]] ){
                        int swap = PQ[front + 2*pqIndex+1];
                        PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                        PQ[front + pqIndex] = swap;
                        pqIndex = 2*pqIndex+1;
                    }
                    else if(Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+2]] && Cx[PQ[front+2*pqIndex+2]] <= Cx[PQ[front+2*pqIndex+1]] ){
                        int swap = PQ[front + 2*pqIndex+2];
                        PQ[front + 2*pqIndex+2] = PQ[front +pqIndex];
                        PQ[front + pqIndex] = swap;
                        pqIndex = 2*pqIndex+2;                    
                    } 
                    else{
                        break;
                    }
                }
                else if(ln[PQ[front+pqIndex]] == ln[PQ[front+2*pqIndex+1]] && Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+1]]){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else if(ln[PQ[front+pqIndex]] == ln[PQ[front+2*pqIndex+2]] && Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+2]]){
                    int swap = PQ[front + 2*pqIndex+2];
                    PQ[front + 2*pqIndex+2] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+2; 
                }
                else{
                    break;
                }
            }
            
        }

        //removed from openList
        openList[node] = -1;

        //added to expand next
        int len = atomicAdd(expandNodes_size,1);
        expandNodes[len]=node;
    }

} 



int main(){

    //the K PQ
    int K ;
    scanf("%d\n",&K);

    FILE* fgraph = fopen("graph_cg.txt","r");
    int N,E;
   
    fscanf(fgraph,"%d %d\n",&N,&E);

    int endNode = N-1;

    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*E);
    //weight is hop count =1
    
    float* H_hx = (float*)malloc(sizeof(float)*N);
    float* H_cx = (float*)malloc(sizeof(float)*N);
    float* H_ln = (float*)malloc(sizeof(float)*N);

    int* H_parent = (int*)malloc(sizeof(int)*N);
    int* H_route = (int*)malloc(sizeof(int)*N);

    int* H_openList = (int*)malloc(sizeof(int)*N);
    int* H_PQ_size = (int*)malloc(sizeof(int)*K);


    //to compute distance
    float* H_cord_x = (float*)malloc(sizeof(float)*N);
    float* H_cord_y = (float*)malloc(sizeof(float)*N);

    float* H_res_energy = (float*)malloc(sizeof(float)*N);

    memset(H_PQ_size,0,sizeof(int)*K);
    memset(H_openList,-1,sizeof(int)*N);
    memset(H_ln,0,sizeof(float)*N);
   
    memset(H_res_energy,E_INIT,sizeof(float)*N);
    
    memset(H_cx,FLT_MAX,sizeof(float)*N);
    
    memset(H_parent,-1,sizeof(int)*N);
    memset(H_route,-1,sizeof(int)*N);

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%d",&H_edges[i]);  
    }

    for(int i=0;i<N;i++){
        fscanf(fgraph,"%d",&H_offset[i]);
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%u",&H_weight[i]);
    }

    fclose(fgraph);

    FILE* f_cord = fopen("Cord.txt","r");
    for(int i=0;i<N;i++){
        int x,y; 
        fscanf(fhx,"%f %f\n",&x,&y);
        H_cord_x[i]= x;
        H_cord_y[i] = y;
    }
    fclose(f_cord);

    float dx,dy;
    dx = H_cord_x[endNode];
    dy = H_cord_y[endNode];
    for(int i=0;i<N;i++){
        float x1,y1;
        x1=H_cord_x[i];
        y1= H_cord_y[i]
        H_hx[i] = sqrt((x1-dx)*(x1-dx)+(y1-dy)*(y1-dy));
    }

    printf("[INFO] completed taking input\n");

    //init Host var
    int* H_flagEnd = (int*)malloc(sizeof(int));
    int* H_flagfound = (int*)malloc(sizeof(int));
    int* H_a0 = (int*)malloc(sizeof(int));

    int* H_nV_size = (int*)malloc(sizeof(int));
    int* H_nV = (int*)malloc(sizeof(int)*N);

    //required coz if many tries to add same in diff threads high low lower
    int* H_nVFlag = (int*)malloc(sizeof(int)*N);
    memset(H_nVFlag,-1,sizeof(int)*N);    

    *H_flagEnd = 0;
    *H_flagfound = 0;
    *H_a0 = 0;

    int* H_end_A_star = (int*)malloc(sizeof(int));
    *H_end_A_star = 0;

    //device var
    //graph struture
    int* D_offset;
    int* D_edges ;
    float* D_hx;
    float* D_Cx;
    float* D_ln;
    
    int* D_parent;

    //Priority queue size
    int* D_PQ_size;

    //flag if in openList(contains which PQ)
    int* D_openList;
    //lock for nodes
    int* D_lock;
    //next nodes flag
    int* D_nVFlag;
    //next nodes array to insert PQ
    int* D_nV;
    int* D_nV_size;
    
    //nodes to be expanded ( extracted from PQ )
    int* D_expandNodes;
    int* D_expandNodes_size;
    
    //flag to end while loop and found the destination
    int* D_flagEnd;
    int* D_flagfound;
 
    float* D_cord_x;
    float* D_cord_y;

    float* D_res_energy;

    gpuErrchk ( cudaMalloc(&D_offset,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_edges,sizeof(int)*E) );

    gpuErrchk ( cudaMalloc(&D_hx,sizeof(float)*N) );
    gpuErrchk ( cudaMalloc(&D_Cx,sizeof(float)*N) );
    gpuErrchk ( cudaMalloc(&D_ln,sizeof(float)*N) );

    gpuErrchk ( cudaMalloc(&D_parent,sizeof(int)*N) );
    
    gpuErrchk ( cudaMalloc(&D_PQ_size,sizeof(int)*K) );
    
    gpuErrchk ( cudaMalloc(&D_openList,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_lock,sizeof(int)*N) );

    //for next set of vertices to add in PQ
    gpuErrchk ( cudaMalloc(&D_nV,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_nV_size,sizeof(int)) );
    gpuErrchk ( cudaMalloc(&D_nVFlag,sizeof(int)*N) );

    //next nodes to expand
    gpuErrchk ( cudaMalloc(&D_expandNodes,sizeof(int)*K) );  //changed to K
    gpuErrchk ( cudaMalloc(&D_expandNodes_size,sizeof(int)) );

    //flag to end search
    gpuErrchk( cudaMalloc(&D_flagEnd,sizeof(int)) );
    gpuErrchk( cudaMalloc(&D_flagfound,sizeof(int)) );
    
    //energy and cords
    gpuErrchk ( cudaMalloc(&D_res_energy,sizeof(float)*N) );
    gpuErrchk ( cudaMalloc(&D_cord_x,sizeof(float)*N) );
    gpuErrchk ( cudaMalloc(&D_cord_y,sizeof(float)*N) );

    //copy
    gpuErrchk ( cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_edges,H_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_hx,H_hx,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_Cx,H_cx,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_ln,H_ln,sizeof(float)*N,cudaMemcpyHostToDevice) );
    
    gpuErrchk ( cudaMemset(D_lock,0,sizeof(int)*N) );

    gpuErrchk ( cudaMemcpy(D_res_energy,H_res_energy,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_cord_x,H_cord_x,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_cord_y,H_cord_y,sizeof(float)*N,cudaMemcpyHostToDevice) );



    int count_round = 0;

    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int N_numBlocks = (N+numThreads-1)/numThreads;

    while(*H_end_A_star == 0){
        //init sol array
        memset(H_route,-1,sizeof(int)*N);

        for(int i=0;i<N;i++){
            //A* algo
            int startNode = i;
            H_cx[startNode]=H_hx[startNode];
            // H_PQ[0]=startNode;
            H_PQ_size[0]=1;
            H_openList[startNode]=0;

            //reset
            gpuErrchk ( cudaMemcpy(D_Cx,H_cx,sizeof(float)*N,cudaMemcpyHostToDevice) );          
            gpuErrchk ( cudaMemcpyToSymbol(PQ,&startNode, sizeof(int), 0, cudaMemcpyHostToDevice) );

            gpuErrchk ( cudaMemcpy(D_parent,H_parent,sizeof(int)*N,cudaMemcpyHostToDevice) );
            gpuErrchk ( cudaMemcpy(D_openList,H_openList,sizeof(int)*N,cudaMemcpyHostToDevice) );
            gpuErrchk ( cudaMemcpy(D_PQ_size,H_PQ_size,sizeof(int)*K,cudaMemcpyHostToDevice) );
        
            gpuErrchk ( cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice) );
            gpuErrchk ( cudaMemcpy(D_flagfound,H_flagfound,sizeof(int),cudaMemcpyHostToDevice) );
           
            gpuErrchk ( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );
        
            gpuErrchk ( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
            gpuErrchk ( cudaMemcpy(D_expandNodes_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
           
            int flag_PQ_not_empty = 0;
            for(int i=0;i<K;i++){
                if(H_PQ_size[i]>0)
                    flag_PQ_not_empty=1;
            }

             //DO A* initailly on whole graph
            while(*H_flagEnd==0 && flag_PQ_not_empty==1){
                
                //extract min
                extractMin<<<numBlocks,numThreads>>>(D_PQ_size, D_expandNodes,D_expandNodes_size,D_Cx,D_openList,N,K);
                
                gpuErrchk(cudaPeekAtLastError() );

                cudaDeviceSynchronize();

                
                A_star_expand<<<numBlocks,numThreads>>>(D_offset,D_edges,D_weight,D_hx,D_parent,D_Cx,
                    D_expandNodes,D_expandNodes_size, D_lock ,D_flagfound,D_openList,
                    N,E,K,endNode,D_nVFlag,D_PQ_size,
                    false,D_diff_offset,D_diff_edges,D_diff_weight,0);
                
                gpuErrchk(cudaPeekAtLastError() );
                
                cudaDeviceSynchronize();


                keepHeapPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_Cx,N,K);
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
                
                //gen from flag D_nV
                //for N in parallel
                setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
                
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
                

                insertPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_nV,D_nV_size,D_Cx,K,N,D_openList);
                
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
            
                //cpy flagend and flagEmpty
                gpuErrchk( cudaMemcpy(H_flagfound,D_flagfound, sizeof(int),cudaMemcpyDeviceToHost) );
                gpuErrchk( cudaMemcpy(H_PQ_size,D_PQ_size, sizeof(int)*K,cudaMemcpyDeviceToHost) );
                
                //reset nVFlag
                gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

                //reset next insert array
                gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
                gpuErrchk( cudaMemcpy(D_expandNodes_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
                

                flag_PQ_not_empty = 0;
                for(int i=0;i<K;i++){
                    if(H_PQ_size[i]>0)
                        flag_PQ_not_empty=1;
                }

                //check for mins
                if( *H_flagfound==1 && flag_PQ_not_empty==1){
                    //end 
                    gpuErrchk( cudaMemcpy(D_flagEnd,H_flagfound,sizeof(int),cudaMemcpyHostToDevice) );

                    checkMIN<<< numBlocks,numThreads >>>(D_PQ_size,D_flagEnd,D_Cx,endNode,N,K);
                    
                    gpuErrchk(cudaPeekAtLastError() );
                    cudaDeviceSynchronize();
                    gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );
                }

            }


        }

        //broadcasted sol

        // update energy l(n)


        count_round++;
    }

    printf("rounds: %d\n",count_round);

    return 0;
}