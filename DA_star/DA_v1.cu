#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  MAX_NODE  1000000
#define  DEBUG 1

__device__ volatile int Cx[MAX_NODE];
__device__ volatile int PQ[MAX_NODE];


//K in parallel
__global__ void extractMin(int* PQ_size, int* expandNodes,int* expandNodes_size,int* openList,int N,int K){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id<K && PQ_size[id]>0){

        //extract min from PQ
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];

        // printf("extract min %d %d\n",id,node);

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
                    pqIndex = 2*pqIndex+1;                    
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
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagEnd,int* openList,
    int N,int E, int K,int dest,int* nVFlag,int* PQ_size,
    int flagDiff,int* diff_off,int* diff_edge,int* diff_weight ){
       
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id< *expandNodes_size ){

        int node = expandNodes[id];
        printf("%d %d\n",id,node);
        //reach dest
        if(node == dest){
            *flagEnd = 1;
            printf("found %d\n",id);
            return;
        }

        // expand
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
            
        while(start<end){ 
            int child = edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            printf("%d$ before while\n",id);

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            while(atomicCAS(&lock[child],0,1)!=0){
            }

            printf("%d$%d: %d ,%d\n",node,child,Cx[child],lock[child]);
            //update cost value
            if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                __threadfence();
                parent[child] = node;

                if(DEBUG)
                    printf("exp: %d %d\n",node,child);

                if(openList[child]>=0){
                    //update operating on one thread
                    if(DEBUG)
                        printf("upd: %d %d\n",node,child);
                    
                    int Kq = openList[child];
                    int front = Kq*( (N+K-1)/K );
                    int index = -1;
                    for(int i=front;i<front+PQ_size[Kq];i++){
                        if(PQ[i]==child){
                            index = i;
                        }
                    }

                    if(index > 0){
                        int i = index;
                        while(i > front){
                            if( Cx[PQ[(i-1)/2]] > Cx[PQ[i]] ){
                                int swap = PQ[i];
                                PQ[i] = PQ[(i-1)/2];
                                PQ[(i-1)/2] = swap;
                                i = (i-1)/2;
                            }
                            else
                                break;
                        }
                    }

                    __threadfence();

                }else{
                    nVFlag[child]=1;
                    //add only once
                }

                
                
            }
            //unlock       
            atomicCAS(&lock[child],1,0);    
            
            start++;
        }

        if(DEBUG)
            printf("%d outside while\n",id);

        if(flagDiff){

            //
            printf("something\n");

            start = diff_off[node];
            end = E;
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
                //loop till acquire the lock
                while(atomicCAS(&lock[child],0,1)!=0){
                }
               
                printf("%d$%d: %d ,%d\n",node,child,Cx[child],lock[child]);
                
                //update cost value
                if( Cx[child] > (Cx[node] - Hx[node])+ diff_weight[start]+ Hx[child] ){
                    Cx[child]  = (Cx[node] - Hx[node])+ diff_weight[start] + Hx[child];
                    __threadfence();
                    parent[child] = node;
                    
                    // printf("%d-%d: %d ,%d\n",node,child,Cx[child],lock[child]);
                    if(openList[child]>=0){
                         //update operating on one thread
                        int Kq = openList[child];
                        int front = Kq*( (N+K-1)/K );
                        int index = -1;
                        for(int i=front;i<front+PQ_size[Kq];i++){
                            if(PQ[i]==child){
                                index = i;
                            }
                        }

                        if(index > 0){
                            int i = index;
                            while(i > front){
                                if( Cx[PQ[(i-1)/2]] > Cx[PQ[i]] ){
                                    int swap = PQ[i];
                                    PQ[i] = PQ[(i-1)/2];
                                    PQ[(i-1)/2] = swap;
                                    i = (i-1)/2;
                                }
                                else
                                    break;
                            }
                        }

                        __threadfence();

                    }else{
                        nVFlag[child]=1;
                    }
                    
                }
                //unlock       
                atomicCAS(&lock[child],1,0);    
                
                start++;
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
        // printf("s: %d %d\n",*nVsize,PQS[id]);
        while(i<*nVsize){
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


__global__ void printCX(int dest){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id==0){
        printf("cost: %d\n",Cx[dest]);
    }
}

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
    int* H_diff_weight = (int*)malloc(sizeof(int)*E);
    
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
        fscanf(fhx,"%d",&H_hx[i]);
    }

    fclose(fgraph);
    fclose(fhx);
    printf("completed input\n");

    //init Host var

    int* H_flagEnd = (int*)malloc(sizeof(int));
    int* H_a0 = (int*)malloc(sizeof(int));

    //required coz if many tries to add same in diff threads high low lower
    int* H_nVFlag = (int*)malloc(sizeof(int)*N);
    memset(H_nVFlag,-1,sizeof(int)*N);    

    *H_flagEnd = 0;
    *H_a0 = 0;

    //insert startNode in PQ[0]
    H_cx[startNode]=H_hx[startNode];
    H_PQ[0]=startNode;
    H_PQ_size[0]=1;
    H_openList[startNode]=0;

    int* D_offset;
    int* D_edges ;
    unsigned int* D_weight;
    int* D_hx;
    int* D_parent;
    // int* D_PQ;
    int* D_PQ_size;

    int* D_openList;
    int* D_lock;

    int* D_diff_edges;
    int* D_diff_offset;
    int* D_diff_weight;

    int* D_nVFlag;
    int* D_nV;
    int* D_nV_size;
    int* D_expandNodes;
    int* D_expandNodes_size;
    int* D_flagEnd;

    cudaMalloc(&D_offset,sizeof(int)*N);
    cudaMalloc(&D_edges,sizeof(int)*E);
    cudaMalloc(&D_weight,sizeof(unsigned int)*E);
    cudaMalloc(&D_hx,sizeof(int)*N);
    cudaMalloc(&D_parent,sizeof(int)*N);
   // cudaMalloc(&D_PQ,sizeof(int)*N);
    cudaMalloc(&D_PQ_size,sizeof(int)*K);
    cudaMalloc(&D_openList,sizeof(int)*N);

    cudaMalloc(&D_lock,sizeof(int)*N);


    //diff csr
    cudaMalloc(&D_diff_edges,sizeof(int)*E);
    cudaMalloc(&D_diff_offset,sizeof(int)*N);
    cudaMalloc(&D_diff_weight,sizeof(int)*E);

    //for next set of vertices to add in PQ
    cudaMalloc(&D_nV,sizeof(int)*N);
    cudaMalloc(&D_nV_size,sizeof(int));
    cudaMalloc(&D_nVFlag,sizeof(int)*N);

    cudaMalloc(&D_expandNodes,sizeof(int)*N);
    cudaMalloc(&D_expandNodes_size,sizeof(int));

    //flag to end search
    cudaMalloc(&D_flagEnd,sizeof(int));
    
    
    cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_edges,H_edges,sizeof(int)*E,cudaMemcpyHostToDevice);
    cudaMemcpy(D_weight,H_weight,sizeof(unsigned int)*E,cudaMemcpyHostToDevice);
    cudaMemcpy(D_hx,H_hx,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_parent,H_parent,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_openList,H_openList,sizeof(int)*N,cudaMemcpyHostToDevice);
    
    cudaMemcpy(D_diff_edges,H_diff_edges,sizeof(int)*E,cudaMemcpyHostToDevice);
    cudaMemcpy(D_diff_offset,H_diff_offset,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_diff_weight,H_diff_weight,sizeof(int)*E,cudaMemcpyHostToDevice);


    // cudaMemcpy(D_PQ,H_PQ,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_PQ_size,H_PQ_size,sizeof(int)*K,cudaMemcpyHostToDevice);


    cudaMemcpyToSymbol(Cx,H_cx, sizeof(int)*N, 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(PQ,H_PQ, sizeof(int)*N, 0, cudaMemcpyHostToDevice);

    cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice);
   
    cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice);

    cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(D_expandNodes_size,H_a0,sizeof(int),cudaMemcpyHostToDevice);

    
    cudaMemset(D_lock,0,sizeof(int)*N);


    int flag_PQ_empty = 0;
    for(int i=0;i<K;i++){
        if(H_PQ_size[i]>0)
            flag_PQ_empty=1;
    }


    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int N_numBlocks = (N+numThreads-1)/numThreads;

    //DO A* initailly on whole graph
    while(*H_flagEnd==0 && flag_PQ_empty==1){
        
        
        //extract min
        extractMin<<<numBlocks,numThreads>>>(D_PQ_size, D_expandNodes,D_expandNodes_size,D_openList,N,K);
        
        cudaDeviceSynchronize();
        if(DEBUG)
            printf("extract min complete\n");


        A_star_expand<<<numBlocks,numThreads>>>(D_offset,D_edges,D_weight,D_hx,D_parent,
            D_expandNodes,D_expandNodes_size, D_lock ,D_flagEnd,D_openList,
            N,E,K,endNode,D_nVFlag,D_PQ_size,
            false,D_diff_offset,D_diff_edges,D_diff_offset );
        
        cudaDeviceSynchronize();
        if(DEBUG)
            printf("expand complete\n");
        //gen from flag D_nV
        //for N in parallel
        setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);

        cudaDeviceSynchronize();
        if(DEBUG)
            printf("set complete\n");

        insertPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_nV,D_nV_size,K,N,D_openList);
        
        cudaDeviceSynchronize();
        if(DEBUG)
            printf("insert complete\n");

        //cpy flagend and flagEmpty
        cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(H_PQ_size,D_PQ_size, sizeof(int)*K,cudaMemcpyDeviceToHost);
        
        //reset nVFlag
        cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice);

        //reset next insert array
        cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice);
        cudaMemcpy(D_expandNodes_size,H_a0,sizeof(int),cudaMemcpyHostToDevice);

        if(DEBUG)
            printf("to host copy complete\n");
        
        flag_PQ_empty = 0;
        for(int i=0;i<K;i++){
            if(H_PQ_size[i]>0)
                flag_PQ_empty=1;
        }
   
    }

    printCX<<<1,1>>>(endNode);
    cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost);


    if(*H_flagEnd==1){
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
