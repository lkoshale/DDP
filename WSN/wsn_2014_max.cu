#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>




/***all macros**/
#define E_INIT  5            // in joules
#define E_ELEC  50e-9          //in nj = 1e-9j
#define E_AMP   100e-12          // in pj = 1e-12j
#define ALPHA   0.7     
#define BETA    0.3
#define DELTA   0.0 


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

//all fn declarations here
void calc_hx(int* offset,int* edges,float* hx,float* cord_x,float* cord_y,int N,int E,int sink);

void update_energy(int* route,float* res_energy,int* Nt,int* Nr,float* cord_x,float* cord_y,int N,int sink);

void check_failure(float* res_energy,int N, int* flag);

/**** device Code *******/

// __device__ volatile int Cx[MAX_NODE];
__device__ volatile int PQ[MAX_NODE];


//K in parallel
__global__ void extractMax(int* PQ_size, int* expandNodes,int* expandNodes_size,float* Cx,int* openList,int N,int K){
    
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
                if( Cx[PQ[front+pqIndex]] < Cx[PQ[front+2*pqIndex+1]]){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else
                    break;
            }
            else{
                if( Cx[PQ[front+pqIndex]] < Cx[PQ[front+2*pqIndex+1]] && Cx[PQ[front+2*pqIndex+1]] >= Cx[PQ[front+2*pqIndex+2]] ){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else if(Cx[PQ[front+pqIndex]] < Cx[PQ[front+2*pqIndex+2]] && Cx[PQ[front+2*pqIndex+2]] >= Cx[PQ[front+2*pqIndex+1]] ){
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

//for K in parallel
__global__ void A_star_expand(int* off,int* edge,float* Hx,int* parent,volatile float* Cx,
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagfound,int* openList,
    int N,int E, int K,int dest,int* nVFlag,int* PQ_size,
    float* res_energy,int* Nt,int* Nr){
       
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id< *expandNodes_size ){

        int node = expandNodes[id];
      
        
        //reach dest
        if(node == dest){
            atomicOr(flagfound,1);
        }

        // expand
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start < end){ 
            int child = edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            bool leaveLoop = false;

            while(leaveLoop==false){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section
                    float exp_weight = ( ALPHA*(res_energy[child]/E_INIT) + BETA*(Nr[child]/Nt[child]) + DELTA*0) ;
                   
                    if( Cx[child] < (Cx[node] - Hx[node])+ exp_weight + Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ exp_weight + Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        if(openList[child]==-1){
                            nVFlag[child]=1;
                            //add only once
                        }
                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);

                }

                __syncthreads();

            }

            start++;
        }
        
    }//end 

}

//K in parallel -- O(N)
__global__ void keepHeapPQ(int* PQ_size,float* Cx,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K && PQ_size[id] > 0){
        int front  = id*( (N+K-1)/K );
        int size = PQ_size[id];
        
        for(int i=front;i<front+size;i++){
            if(2*i+2 < front+size){
                int cost = Cx[PQ[i]];
                int costLeft = Cx[PQ[2*i+1]];
                int costRight = Cx[PQ[2*i+2]]; 

                if( cost < costLeft  ||  cost < costRight  ){
                    int index ;
                    if(costLeft >= costRight)
                        index = 2*i+1;
                    else
                        index = 2*i+2;
                    
                    while(index > front){
                        if( Cx[PQ[(index-1)/2]] < Cx[PQ[index]] ){
                            int swap = PQ[index];
                            PQ[index] = PQ[(index-1)/2];
                            PQ[(index-1)/2] = swap;
                            index = (index-1)/2;
                        }
                        else
                            break;
                    }
                }
            }
            else if(2*i+1 < front+size){
                if(Cx[PQ[i]] < Cx[PQ[2*i+1]]){
                    int index = 2*i+1;
                    while(index > front){
                        if( Cx[PQ[(index-1)/2]] < Cx[PQ[index]] ){
                            int swap = PQ[index];
                            PQ[index] = PQ[(index-1)/2];
                            PQ[(index-1)/2] = swap;
                            index = (index-1)/2;
                        }
                        else
                            break;
                    }
                }
            }
        }
    }
}

//N threads
__global__ void setNV(int* nextFlag,int* nextV,int* nvSize,int N){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < N){
        if(nextFlag[id]==1){
            int index = atomicAdd(nvSize,1);
            nextV[index]=id;
        }
    }
}


//for K in parallel
__global__ void insertPQ(int* PQS,int* nextV,int* nVsize,float* Cx,int K,int N,int* openList){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K){

        int front = id*( (N+K-1)/K );
        int i = id;
        
        while(i<*nVsize){            
            //if not already present
            if(openList[nextV[i]]!=-1){
                i+=K;
                continue;
            }

            PQ[front+PQS[id]]= nextV[i];
            PQS[id]+=1;

            //add in openList
            openList[nextV[i]] = id;

            if(PQS[id]>1){
                int index = PQS[id]-1;
                while(index>0){
                    if(Cx[PQ[front+ (index-1)/2]] < Cx[PQ[front+index]]){
                        int swap = PQ[front+index];
                        PQ[front+index]=PQ[front+ (index-1)/2];
                        PQ[front+ (index-1)/2] = swap;
                        index = (index-1)/2;
                    }
                    else
                        break;
                }
            }
            i += K;
        }
    }
}


//for K in parallel
__global__ void checkMAX(int* PQ_size,int* flagEnd,float* Cx,int dest,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < K && PQ_size[id] > 0 ){
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];
        //check if atleast one max, dont end the a*
        if( Cx[node] > Cx[dest] ){
            atomicAnd(flagEnd,0);
        }
    }
}


int main(){

    //the K PQ
    int K ;
    scanf("%d",&K);

    FILE* fgraph = fopen("graph.txt","r");
    int N,E;
   
    fscanf(fgraph,"%d %d\n",&N,&E);

    int endNode = N-1;

    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*E);
    //weight is hop count =1
    
    float* H_hx = (float*)malloc(sizeof(float)*N);
    float* H_cx = (float*)malloc(sizeof(float)*N);

    int* H_parent = (int*)malloc(sizeof(int)*N);
    int* H_route = (int*)malloc(sizeof(int)*N);

    int* H_openList = (int*)malloc(sizeof(int)*N);
    int* H_PQ_size = (int*)malloc(sizeof(int)*K);


    //to compute distance
    float* H_cord_x = (float*)malloc(sizeof(float)*N);
    float* H_cord_y = (float*)malloc(sizeof(float)*N);

    float* H_res_energy = (float*)malloc(sizeof(float)*N);

    int* H_packet_recv =  (int*)malloc(sizeof(int)*N);
    int* H_packet_sent = (int*)malloc(sizeof(int)*N);

    memset(H_PQ_size,0,sizeof(int)*K);
    memset(H_openList,-1,sizeof(int)*N);
    
    memset(H_parent,-1,sizeof(int)*N);
    memset(H_route,-1,sizeof(int)*N);

    for(int i=0;i<N;i++){
        H_cx[i] = 0;
        H_res_energy[i] = E_INIT;
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%d",&H_edges[i]);  
    }

    for(int i=0;i<N;i++){
        fscanf(fgraph,"%d",&H_offset[i]);
    }

    fclose(fgraph);

    FILE* f_cord = fopen("Cord.txt","r");
    for(int i=0;i<N;i++){
        float x,y; 
        fscanf(f_cord,"%f %f\n",&x,&y);
        H_cord_x[i]= x;
        H_cord_y[i] = y;
    }
    fclose(f_cord);

    calc_hx(H_offset,H_edges,H_hx,H_cord_x,H_cord_y,N,E,endNode);

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
    int* D_packet_recv;
    int* D_packet_sent;

    gpuErrchk ( cudaMalloc(&D_offset,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_edges,sizeof(int)*E) );

    gpuErrchk ( cudaMalloc(&D_hx,sizeof(float)*N) );
    gpuErrchk ( cudaMalloc(&D_Cx,sizeof(float)*N) );
   
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
    gpuErrchk( cudaMalloc(&D_packet_recv,sizeof(int)*N) );
    gpuErrchk( cudaMalloc(&D_packet_sent,sizeof(int)*N) );

    //copy
    gpuErrchk ( cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_edges,H_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_hx,H_hx,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_Cx,H_cx,sizeof(float)*N,cudaMemcpyHostToDevice) );
    
    gpuErrchk ( cudaMemset(D_lock,0,sizeof(int)*N) );

    gpuErrchk ( cudaMemcpy(D_res_energy,H_res_energy,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_cord_x,H_cord_x,sizeof(float)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_cord_y,H_cord_y,sizeof(float)*N,cudaMemcpyHostToDevice) );
    
    //1 TO EVADE DIVIDE BY 0 ERROR
    gpuErrchk ( cudaMemset(D_packet_recv,1,sizeof(int)*N) );
    gpuErrchk ( cudaMemset(D_packet_sent,1,sizeof(int)*N) );


    int count_round = 0;

    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int N_numBlocks = (N+numThreads-1)/numThreads;

    while(*H_end_A_star == 0){
        //init sol array
        memset(H_route,-1,sizeof(int)*N);

        for(int i=0;i<N;i++){

            if(H_route[i]!=-1 || i == endNode)
                continue;
            
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
                extractMax<<<numBlocks,numThreads>>>(D_PQ_size, D_expandNodes,D_expandNodes_size,D_Cx,D_openList,N,K);
                
                gpuErrchk(cudaPeekAtLastError() );

                cudaDeviceSynchronize();
                
                A_star_expand<<<numBlocks,numThreads>>>(D_offset,D_edges,D_hx,D_parent,D_Cx,
                    D_expandNodes,D_expandNodes_size, D_lock ,D_flagfound,D_openList,
                    N,E,K,endNode,D_nVFlag,D_PQ_size,
                    D_res_energy,D_packet_sent,D_packet_recv);
                
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

                    checkMAX<<< numBlocks,numThreads >>>(D_PQ_size,D_flagEnd,D_Cx,endNode,N,K);
                    
                    gpuErrchk(cudaPeekAtLastError() );
                    cudaDeviceSynchronize();
                    gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );
                }

            }

            gpuErrchk( cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );
            int p = endNode;
            printf("PATH: %d ",endNode);
            while(H_parent[p]!=-1){
                
                //create route
                H_route[H_parent[p]] = p;

                p = H_parent[p];
                printf("%d ",p);
                if(p==startNode)
                    break;
            }
            printf("\n");


            memset(H_PQ_size,0,sizeof(int)*K);
            memset(H_openList,-1,sizeof(int)*N);
            memset(H_cx,0,sizeof(float)*N);
            memset(H_parent,-1,sizeof(int)*N);

            *H_flagEnd = 0;
            *H_flagfound = 0;
            
        }

        //broadcasted sol

        // update energy
        update_energy(H_route,H_res_energy,H_packet_sent,H_packet_recv,H_cord_x,H_cord_y,N,endNode);        

        gpuErrchk ( cudaMemcpy(D_packet_recv,H_packet_recv,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_packet_sent,H_packet_sent,sizeof(int)*N,cudaMemcpyHostToDevice) );

        gpuErrchk ( cudaMemcpy(D_res_energy,H_res_energy,sizeof(float)*N,cudaMemcpyHostToDevice) );

        //check for end
        check_failure(H_res_energy,N,H_end_A_star);

        count_round++;
    }

    printf("rounds: %d\n",count_round);

    return 0;
}

void calc_hx(int* offset,int* edges,float* hx,float* cord_x,float* cord_y,int N,int E,int sink){
    float dx,dy;
    dx = cord_x[sink];
    dy = cord_y[sink];

    for(int i=0;i<N;i++){
        float x,y;
        x = cord_x[i];
        y = cord_y[i];    

        int start = offset[i];
        int end = N;
        if(i!=N-1)
            end = offset[i+1];
        
        float sum = 0;
        int count = 0;

        while(start < end){
            int child = edges[start];
            sum+= sqrtf( (x-cord_x[child])*(x-cord_x[child]) + (y-cord_y[child])*(y-cord_y[child]) );
            start++;
            count++;
        }

        float dist = sqrtf( (x-dx)*(x-dx) + (y-dy)*(y-dy) );
        //hop counts 
        if(dist!=0){
            hx[i] =  sum / (dist * count ) ;
        }
        else
            hx[i] = 0;
    }

}

void update_energy(int* route,float* res_energy,int* Nt,int* Nr,float* cord_x,float* cord_y,int N,int sink){

    int k = 2048;    //256 bit packet length

    for(int j=0;j<N;j++){
        int start = j;
        while(route[start]!=-1){
            int node = route[start];
            Nr[node]+=1;
            res_energy[node]-= k*E_ELEC;
            
            Nt[start]+=1;
            float dist_sq = (cord_x[start]-cord_x[node])*(cord_x[start]-cord_x[node]) +(cord_y[start]-cord_y[node])*(cord_y[start]-cord_y[node]);
            float energy = k*(E_ELEC + E_AMP*dist_sq);
            res_energy[start]-=energy*100;
           
            start = node;
            
        }

        if(start != sink)
            printf("ERROR; invalid path\n");

    }

}

void check_failure(float* res_energy,int N, int* flag){
    for(int i=0;i<N;i++){
       // printf("%d:%f\n",i,res_energy[i]);
        if(res_energy[i]<0){
            *flag =1;
            printf("dead: %d\n",i);
            break;
        }       
    }
}