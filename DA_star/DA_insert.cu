/*
 For DIRECTED GRAPH
*/

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>


/***all macros**/
#define  MAX_NODE  100000000
#define  DEBUG 0

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

/***function declarations***/
void insertDiff(unordered_map< unsigned int, Node*>& Graph,int a,int b,unsigned int c);
void createDiffGraph(int N,unordered_map<unsigned int,Node*>& Graph,
    int* diffOff,int* diffEdges,unsigned int* diffWeight );


/**** device Code *******/

__device__ volatile int Cx[MAX_NODE];
__device__ volatile int PQ[MAX_NODE];


//K in parallel
__global__ void extractMin(int* PQ_size, int* expandNodes,int* expandNodes_size,int* openList,int N,int K){
    
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
                if( Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+1]]){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else
                    break;
            }
            else{
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
            
        }

        //removed from openList
        openList[node] = -1;

        //added to expand next
        int len = atomicAdd(expandNodes_size,1);
        expandNodes[len]=node;
    }

} 


//for K in parallel
__global__ void A_star_expand(int* off,int* edge,unsigned int* W,int* Hx,int* parent,
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagfound,int* openList,
    int N,int E, int K,int dest,int* nVFlag,int* PQ_size,
    int flagDiff,int* diff_off,int* diff_edge,int* diff_weight,int dE ){
       
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id< *expandNodes_size ){

        int node = expandNodes[id];
        
        //reach dest
        if(node == dest){
            atomicOr(flagfound,1);
           // *flagfound = 1;
            printf("found %d\n",id);
        }

        // expand
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        
        if(DEBUG)
            printf("%d$expand %d:\n",id,start);

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
                    if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        if(DEBUG)
                            printf("exp: %d %d\n",node,child);
        
                        if(openList[child]==-1){
                            nVFlag[child]=1;
                            //add only once
                        }
                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);

                    // if(DEBUG)
                    //     printf("im not deadlocked\n");
                }

                __syncthreads();
                if(DEBUG)
                    printf("%d$ stuck here\n",id);
                
            }

            start++;
        }

        //diff expand
        if(flagDiff){

            start = diff_off[node];
            end = dE;
            if(node!=N-1)
                end = diff_off[node+1];

            while(start<end){ 
                int child = diff_edge[start];
                
                //deleted edges
                if(child<0){
                    start++;
                    continue;
                }
    
                //array L initilaized with 0
                //get the lock for child to update C(x)
                bool leaveLoop = false;

                while(!leaveLoop){

                    if(atomicCAS(&lock[child],0,1)==0){
                        //critical section
                        if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                            Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                            __threadfence();
                            parent[child] = node;
            
                            if(DEBUG)
                                printf("exp: %d %d\n",node,child);
            
                            if(openList[child]==-1){
                                nVFlag[child]=1;
                                //add only once
                            }
                        }

                        //end critical section
                        leaveLoop = true;

                        atomicCAS(&lock[child],1,0);

                        // if(DEBUG)
                        //     printf("im not deadlocked\n");
                    }

                    __syncthreads();
                    
                }
                
                start++;
            }
            
        }
        //end diff
        
    }//end 

}


//K in parallel -- O(N)
__global__ void keepHeapPQ(int* PQ_size,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K && PQ_size[id] > 0){
        int front  = id*( (N+K-1)/K );
        int size = PQ_size[id];
        
        for(int i=front;i<front+size;i++){
            if(2*i+2 < front+size){
                int cost = Cx[PQ[i]];
                int costLeft = Cx[PQ[2*i+1]];
                int costRight = Cx[PQ[2*i+2]]; 
                if( cost > costLeft  ||  cost > costRight  ){
                    int index ;
                    if(costLeft <= costRight)
                        index = 2*i+1;
                    else
                        index = 2*i+2;
                    
                    while(index > front){
                        if( Cx[PQ[(index-1)/2]] > Cx[PQ[index]] ){
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
                if(Cx[PQ[i]] > Cx[PQ[2*i+1]]){
                    int index = 2*i+1;
                    while(index > front){
                        if( Cx[PQ[(index-1)/2]] > Cx[PQ[index]] ){
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
        //printf("2: %d %d\n",id,nextFlag[id]);
        if(nextFlag[id]==1){
            int index = atomicAdd(nvSize,1);
            nextV[index]=id;
          //  printf("2: %d\n",id);
        }
    }
}


//for K in parallel
__global__ void insertPQ(int* PQS,int* nextV,int* nVsize,int K,int N,int* openList){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K){
        // printf("id: %d\n",id);
        int front = id*( (N+K-1)/K );
        int i = id;
        // if(id==0)
        //     printf("insert: %d\n",*nVsize);
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

            //printf("insert: %d, %d\n",nextV[i],PQS[id]);
            if(PQS[id]>1){
                int index = PQS[id]-1;
                while(index>0){
                    if(Cx[PQ[front+ (index-1)/2]] > Cx[PQ[front+index]]){
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
__global__ void checkMIN(int* PQ_size,int* flagEnd,int dest,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < K && PQ_size[id] > 0 ){
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];
        //check if atleast one min, dont end the a*
       // printf("%d ",Cx[node]);
        if(Cx[dest] > Cx[node] ){
            atomicAnd(flagEnd,0);
        }
    }
}


//add inserted edges to propogate
__global__ void propogateAdd(int* diff_off, int* diff_edges,unsigned int* diff_W,int* Hx,int* addFlag,
            int* lock, int* parent, int N, int dE){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < N){
        int node = id;
        
        int start = diff_off[node];
        int end = dE;
        if(node!=N-1)
            end = diff_off[node+1];
        
        while(start < end ){
            int child = diff_edges[start];
            if(DEBUG)
                printf("ch:%d,%d\n",child,Cx[child]);
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            bool leaveLoop = false;

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section
                    if( Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        
                        addFlag[child]=1;
                        if(DEBUG)
                            printf("add:%d\n",child);

                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);
                }

                __syncthreads();
            }

            start++;
        }
    }    

}

//propogate the change
__global__ void propogate(int* nodes, int* size, int* off, int* edge,unsigned int* W,int* Hx,
                    int N,int E, int* lock, int* parent,int* addFlag,
                    int* diff_off,int* diff_edge,unsigned int* diff_W,int dE){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < *size){
        int node = nodes[id];
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start < end ){
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

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section
                    if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);
                }

                __syncthreads();
                
            }

            start++;
        }    


        start = diff_off[node];
        end = dE;
        if(node!=N-1)
            end = diff_off[node+1];
        
        while(start < end ){
            int child = diff_edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            bool leaveLoop = false;

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section
                    if( Cx[child] != INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);
                }

                __syncthreads();
                
            }

            start++;
        }    

    }
                        
}


//do in 1 thread
__global__ void insertDest(int* PQ_size, int dest,int* openList){
    int id = 0;
    int front = 0;
    if(openList[dest]==-1){
        PQ[front+PQ_size[id]]= dest;
        PQ_size[id]+=1;

        //add in openList
        openList[dest] = id;

        if(PQ_size[id]>1){
            int index = PQ_size[id]-1;
            while(index>0){
                if(Cx[PQ[front+ (index-1)/2]] > Cx[PQ[front+index]]){
                    int swap = PQ[front+index];
                    PQ[front+index]=PQ[front+ (index-1)/2];
                    PQ[front+ (index-1)/2] = swap;
                    index = (index-1)/2;
                }
                else
                    break;
            }
        }
    }
    
}

__global__ void getCx(int dest,int* val){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id==0){
        printf("cost: %d\n",Cx[dest]);
        *val = Cx[dest];
    }
}



/**** main function ****/
int main(){

    //the K PQ
    int K ;
    scanf("%d\n",&K);

    int startNode,endNode;
    scanf("%d %d",&startNode,&endNode);

    FILE* fgraph = fopen("graph.txt","r");
    int N,E;
    fscanf(fgraph,"%d %d\n",&N,&E);

    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*E);
    unsigned int* H_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);
    int* H_hx = (int*)malloc(sizeof(int)*N);
    int* H_cx = (int*)malloc(sizeof(int)*N);
    int* H_parent = (int*)malloc(sizeof(int)*N);


    int* H_PQ = (int*)malloc(sizeof(int)*N);
    int* H_openList = (int*)malloc(sizeof(int)*N);
    int* H_PQ_size = (int*)malloc(sizeof(int)*K);


    //for diff
    int* H_diff_edges = (int*)malloc(sizeof(int)*E);
    int* H_diff_offset = (int*)malloc(sizeof(int)*N);
    unsigned int* H_diff_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);

    //for cost of endNode
    int* H_dest_cost = (int*)malloc(sizeof(int));
    
    memset(H_PQ_size,0,sizeof(int)*K);
    memset(H_parent,-1,sizeof(int)*N);
    memset(H_openList,-1,sizeof(int)*N);

    //init cx
    for(int i=0;i<N;i++){
        H_cx[i]=INT_MAX;
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%d",&H_edges[i]);
    }

    for(int i=0;i<N;i++){
        fscanf(fgraph,"%d",&H_offset[i]);
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%u",&H_weight[i]);
    }

    FILE* fhx = fopen("Hx.txt","r");
    for(int i=0;i<N;i++){
        int temp; 
        fscanf(fhx,"%d",&temp);
        if(temp!=-1)
            H_hx[i]= temp;
        else 
            H_hx[i] = 0;   //to change
    }

    fclose(fgraph);
    fclose(fhx);
    printf("completed input\n");

    //init Host var
    int* H_flagEnd = (int*)malloc(sizeof(int));
    int* H_flagfound = (int*)malloc(sizeof(int));
    int* H_a0 = (int*)malloc(sizeof(int));

    int* H_nV_size = (int*)malloc(sizeof(int));

    //required coz if many tries to add same in diff threads high low lower
    int* H_nVFlag = (int*)malloc(sizeof(int)*N);
    memset(H_nVFlag,-1,sizeof(int)*N);    

    *H_flagEnd = 0;
    *H_flagfound = 0;
    *H_a0 = 0;

    //insert startNode in PQ[0]
    H_cx[startNode]=H_hx[startNode];
    H_PQ[0]=startNode;
    H_PQ_size[0]=1;
    H_openList[startNode]=0;


    //graph struture
    int* D_offset;
    int* D_edges ;
    unsigned int* D_weight;
    int* D_hx;
    int* D_parent;
    
    //Priority queue size
    int* D_PQ_size;

    //flag if in openList(contains which PQ)
    int* D_openList;
    //lock for nodes
    int* D_lock;

    //Diff structure
    int* D_diff_edges;
    int* D_diff_offset;
    unsigned int* D_diff_weight;

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

    //cost of endNode
    int* D_dest_cost;

    gpuErrchk ( cudaMalloc(&D_offset,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_weight,sizeof(unsigned int)*E) );
    gpuErrchk ( cudaMalloc(&D_hx,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_parent,sizeof(int)*N) );
   
    gpuErrchk ( cudaMalloc(&D_PQ_size,sizeof(int)*K) );
    
    gpuErrchk ( cudaMalloc(&D_openList,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_lock,sizeof(int)*N) );

    gpuErrchk ( cudaMalloc(&D_dest_cost,sizeof(int)) );


    //diff csr
    gpuErrchk ( cudaMalloc(&D_diff_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_diff_offset,sizeof(int)*(N+1) ) );   //coz
    gpuErrchk ( cudaMalloc(&D_diff_weight,sizeof(unsigned int)*E) );

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
    
    
    gpuErrchk ( cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_edges,H_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_weight,H_weight,sizeof(unsigned int)*E,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_hx,H_hx,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_parent,H_parent,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_openList,H_openList,sizeof(int)*N,cudaMemcpyHostToDevice) );


    gpuErrchk ( cudaMemcpy(D_PQ_size,H_PQ_size,sizeof(int)*K,cudaMemcpyHostToDevice) );


    gpuErrchk ( cudaMemcpyToSymbol(Cx,H_cx, sizeof(int)*N, 0, cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpyToSymbol(PQ,H_PQ, sizeof(int)*N, 0, cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_flagfound,H_flagfound,sizeof(int),cudaMemcpyHostToDevice) );
   
    gpuErrchk ( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_expandNodes_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );

    
    gpuErrchk ( cudaMemset(D_lock,0,sizeof(int)*N) );


    int flag_PQ_not_empty = 0;
    for(int i=0;i<K;i++){
        if(H_PQ_size[i]>0)
            flag_PQ_not_empty=1;
    }


    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int N_numBlocks = (N+numThreads-1)/numThreads;

    //DO A* initailly on whole graph
    while(*H_flagEnd==0 && flag_PQ_not_empty==1){
        
        //extract min
        extractMin<<<numBlocks,numThreads>>>(D_PQ_size, D_expandNodes,D_expandNodes_size,D_openList,N,K);
        
        gpuErrchk(cudaPeekAtLastError() );

        cudaDeviceSynchronize();

        
        A_star_expand<<<numBlocks,numThreads>>>(D_offset,D_edges,D_weight,D_hx,D_parent,
            D_expandNodes,D_expandNodes_size, D_lock ,D_flagfound,D_openList,
            N,E,K,endNode,D_nVFlag,D_PQ_size,
            false,D_diff_offset,D_diff_edges,D_diff_offset,0);
        
        gpuErrchk(cudaPeekAtLastError() );
        
        cudaDeviceSynchronize();


        keepHeapPQ<<<numBlocks,numThreads>>>(D_PQ_size,N,K);
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        
        //gen from flag D_nV
        //for N in parallel
        setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        

        insertPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_nV,D_nV_size,K,N,D_openList);
        
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

            checkMIN<<< numBlocks,numThreads >>>(D_PQ_size,D_flagEnd,endNode,N,K);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();
            gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );
           // printf("\ninside MIN\n");
        }
   
    }

    getCx<<<1,1>>>(endNode,D_dest_cost);

    if(*H_flagfound==1){
        int p = endNode;
        while(H_parent[p]!=-1){
            printf("%d ",p);
            p = H_parent[p];
        } 
        printf("%d\n",startNode);
    }
    else{
        printf("not found\n");
    }

    

    ///////////////////////////////////////////////
    // A star complete //


    FILE* fdiff = fopen("Updates.txt","r");
    int line;
    while(fscanf(fdiff,"%d\n",&line)!=EOF){
        unordered_map<unsigned int,Node*> Graph;
        int insertEdge=0, delEdge=0;
        for(int i=0;i<line;i++){
            int flag;
            int u,v;
            unsigned int w;
            fscanf(fdiff,"%d %d %d %u\n",&flag,&u,&v,&w);
            if(flag==1){
                insertDiff(Graph,u,v,w);
                insertEdge++;
            }
        }

        createDiffGraph(N,Graph,H_diff_offset,H_diff_edges,H_diff_weight);
        
        if(DEBUG){
            for(int i=0;i<N;i++)
                printf("%d ",H_diff_offset[i]);
            printf("\n");
            for(int i=0;i<insertEdge;i++)
                printf("%d ",H_diff_edges[i]);
            printf("\n");
            
        }

        gpuErrchk ( cudaMemcpy(D_diff_edges,H_diff_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_diff_offset,H_diff_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_diff_weight,H_diff_weight,sizeof(int)*E,cudaMemcpyHostToDevice) );

        //reset D_nV flag
        gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );
    
        //N parallel
        propogateAdd<<<N_numBlocks,numThreads>>>(D_diff_offset, D_diff_edges,D_diff_weight,D_hx,D_nVFlag,
            D_lock,D_parent,N,insertEdge);
        
        //add in next set of nodes
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        

        gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );

        //gen from flag D_nV
        setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();

        // //add in priority queue
        // insertPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_nV,D_nV_size,K,N,D_openList);
        
        // gpuErrchk(cudaPeekAtLastError() );
        // cudaDeviceSynchronize();

        //copy back
        gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );

        //reset nV flags
        gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

        while(*H_nV_size > 0){
            printf("update size: %d\n",*H_nV_size);

            numBlocks = (*H_nV_size+numThreads-1)/numThreads;

            propogate<<<numBlocks,numThreads>>>(D_nV,D_nV_size,D_offset,D_edges,D_weight,D_hx,
                N,E,D_lock,D_parent,D_nVFlag,
                D_diff_offset,D_diff_edges,D_diff_weight,insertEdge);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();

            //reset size=0
            gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );

            //gen from flag D_nV
            setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();

            //copy back
            gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );

            //reset nV flags
            gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );
        }

        //propogate complete do normal A*
        numBlocks = (K+numThreads-1)/numThreads;
        
        //update PQ after propogate
        keepHeapPQ<<<numBlocks,numThreads>>>(D_PQ_size,N,K);
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();

        
        getCx<<<1,1>>>(endNode,D_dest_cost);
        gpuErrchk( cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );


        // found or not found based on Cx
        gpuErrchk( cudaMemcpy(H_dest_cost,D_dest_cost, sizeof(int),cudaMemcpyDeviceToHost) );

        printf("dcost: %d\n",*H_dest_cost);
        if(*H_dest_cost!=INT_MAX){
            int p = endNode;
            while(H_parent[p]!=-1){
                printf("%d ",p);
                p = H_parent[p];
            } 
            printf("%d\n",startNode);
        }
        else{
            printf("not found\n");
        }


    }

}


void insertDiff(unordered_map< unsigned int, Node*>& Graph,int a,int b,unsigned int c){
    unordered_map<unsigned int,Node*>:: iterator itr;
    itr = Graph.find(a);
    if(itr!=Graph.end()){
        Node* n = itr->second;
        unordered_map<unsigned int,Node*>:: iterator it;
        it = Graph.find(b);
        if(it!=Graph.end()){
            Node* v = it->second;
            n->addEdge(v,c);
        }
        else{
            Node* v = new Node(b);
            n->addEdge(v,c);
            Graph.insert(pair<unsigned int,Node*>(b,v));
        }

    }
    else{
        Node* n =new Node(a);
        Graph.insert(pair<unsigned int,Node*>(a,n));

        unordered_map<unsigned int,Node*>:: iterator it;
        it = Graph.find(b);
        if(it!=Graph.end()){
            Node* v = it->second;
            n->addEdge(v,c);
        }
        else{
            Node* v = new Node(b);
            n->addEdge(v,c);
            Graph.insert(pair<unsigned int,Node*>(b,v));
        }

    }
}


void createDiffGraph(int N,unordered_map<unsigned int,Node*>& Graph,
                   int* diffOff,int* diffEdges,unsigned int* diffWeight ){

    int offindex = 0;
    diffOff[offindex] = 0;
    offindex++;
    int k =0;
    int weightCount = 0;

    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->Edges.size();j++){
                diffEdges[k] = n->Edges[j]->val;
                k++;
            }

            for(int j=0;j<n->weights.size();j++){
                diffWeight[weightCount] = n->weights[j];
                weightCount++;
            }

            diffOff[offindex] = k;
            offindex++;
        }
        else{
            diffOff[offindex] = k;
            offindex++;
        }
    }

}
