#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#define  MAX_EDGE  100000000
#define  MAX_NODE  1000000

__device__ volatile int Cx[MAX_NODE];

__global__ void A_star(int* off,int* edge,int* W,int* Hx,int* P,int* PQ,int* PQS,int* L,int* nextFlag,int* endF,int N,int E,int K,int dest){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K && PQS[id]>0){
        //extract min from PQ
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];
        // restructure the heap
        PQ[front]=PQ[front+PQS[id]-1];
        PQS[id]-=1;
        int pqIndex = 0;
        while(2*pqIndex+1 < PQS[id]){
            if(2*pqIndex+2 >= PQS[id]){
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
        //reach dest
        if(node == dest){
            *endF = 1;
        }

        // printf("%d,%d\n",id,PQS[id]);

        // expand
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start<end){
            int child = edge[start];
            // printf("c %d - %d\n",node,child);
            
            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            while(atomicCAS(&L[child],0,1)!=0){
            }
            printf("%d$%d: %d ,%d\n",node,child,Cx[child],L[child]);
            //update cost value
            if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                __threadfence();
                P[child] = node;
                printf("%d-%d: %d ,%d\n",node,child,Cx[child],L[child]);
                nextFlag[child]=1;
            }
            
            //unlock       
            atomicCAS(&L[child],1,0);
            
	          // L[child]=0;

            //printf("cx: %d \n",Cx[child]);
            
            start++;
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
__global__ void insertPQ(int* PQ,int* PQS,int* nextV,int* nVsize,int K,int N){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K){
        // printf("id: %d\n",id);
        int front = id*( (N+K-1)/K );
        int i = id;
        // printf("s: %d %d\n",*nVsize,PQS[id]);
        while(i<*nVsize){
            PQ[front+PQS[id]]= nextV[i];
            PQS[id]+=1;
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

void print_Q(int* H_Kqueue,int* H_Kq_size,int N,int K){
    for(int id=0;id<K;id++){
        printf("%d: ",id);
        int st = id* ( (N+K-1)/K );
        int j = 0;
        while(j < H_Kq_size[id]){
            printf("%d ",H_Kqueue[st+j]);
            j++;
        }
        printf("\n");
    }
}


int main(){
    int N;
    int K;
    scanf("%d %d\n",&N,&K);
    int startNode;
    int endNode;
    scanf("%d %d\n",&startNode,&endNode);


    int* H_hx = (int*)malloc(sizeof(int)*N);
    int* H_cx = (int*)malloc(sizeof(int)*N);

    int* H_parent = (int*)malloc(sizeof(int)*N);

    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*MAX_EDGE);
    int* H_weights = (int*)malloc(sizeof(int)*MAX_EDGE);

    for(int j=0;j<N;j++){
        scanf("%d",&H_hx[j]);
    }

    memset(H_offset,-1,sizeof(int)*N);
    memset(H_cx,100,sizeof(int)*N);
    memset(H_parent,-1,sizeof(int)*N);

    int a,b,c;
    int prev_node = -1;
    int edge_size = 0;
    while(scanf("%d %d %d\n",&a,&b,&c)!=EOF){
        if(a==prev_node){
            H_edges[edge_size]=b;
            H_weights[edge_size]=c;
            edge_size++;
        }
        else{
            H_offset[a]=edge_size;
            H_edges[edge_size]=b;
            H_weights[edge_size]=c;
            edge_size++;
            prev_node = a;

        }
    }

    for(int i=0;i<N;i++){
        if(H_offset[i]==-1){
            int j = i+1;
            int flag = 0;
            while(j<N){
                if(H_offset[j]==-1){
                }
                else{
                    H_offset[i]=H_offset[j];
                    flag= 1;
                    break;
                }
                j++;
            }
            if(flag==0){
                H_offset[i] = edge_size;
            }
        }
    }

    printf("completed input\n");

    // K priority queue with i start at i*(N+K-1)/K to (i+1)(N+K-1)/K or N
    int* H_Kqueue = (int*)malloc(sizeof(int)*N);
    int* H_Kq_size = (int*)malloc(sizeof(int)*K);
    memset(H_Kq_size,0,sizeof(int)*K);
    H_Kqueue[0]=startNode;
    H_cx[startNode]=H_hx[startNode];
    H_Kq_size[0]=1;

    //N lock var
    int* H_lock = (int*)malloc(sizeof(int)*N);
    memset(H_lock,0,sizeof(int)*N);     //unlocked

    int* H_nextFlag = (int*)malloc(sizeof(int)*N);
    memset(H_nextFlag,-1,sizeof(int)*N);     //

    int* H_flagEnd = (int*)malloc(sizeof(int));
    int* a0 = (int*)malloc(sizeof(int));
    *H_flagEnd = 0;     //false
    *a0 = 0;            // to set sizes 0

    int* D_offset;
    int* D_edges;
    int* D_weights;
    int* D_hx;
    int* D_parent;
    int* D_Kqueue;
    int* D_Kq_size;
    int *D_lock;
    int* D_nextFlag;
    int* D_nextV;
    int* D_nV_size;
    int* D_flagEnd;


    cudaMalloc(&D_offset,sizeof(int)*N);
    cudaMalloc(&D_edges,sizeof(int)*edge_size);
    cudaMalloc(&D_weights,sizeof(int)*edge_size);
    cudaMalloc(&D_hx,sizeof(int)*N);
    cudaMalloc(&D_parent,sizeof(int)*N);
    cudaMalloc(&D_Kqueue,sizeof(int)*N);
    cudaMalloc(&D_Kq_size,sizeof(int)*K);
    cudaMalloc(&D_lock,sizeof(int)*N);
    cudaMalloc(&D_nextFlag,sizeof(int)*N);
    cudaMalloc(&D_nextV,sizeof(int)*N);
    cudaMalloc(&D_nV_size,sizeof(int));
    cudaMalloc(&D_flagEnd,sizeof(int));


    cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_edges,H_edges,sizeof(int)*edge_size,cudaMemcpyHostToDevice);
    cudaMemcpy(D_weights,H_weights,sizeof(int)*edge_size,cudaMemcpyHostToDevice);
    cudaMemcpy(D_hx,H_hx,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_parent,H_parent,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_Kqueue,H_Kqueue,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_Kq_size,H_Kq_size,sizeof(int)*K,cudaMemcpyHostToDevice);
    cudaMemcpy(D_lock,H_lock,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_nextFlag,H_nextFlag,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_nV_size,a0,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(Cx,H_cx, sizeof(int)*N, 0, cudaMemcpyHostToDevice);

    int flagQempty = 0;
    for(int i=0;i<K;i++){
        if(H_Kq_size[i]>0)
            flagQempty=1;
    }

    int j = 0;

    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int new_Blocks = (N + numThreads-1)/numThreads;
    
    while(*H_flagEnd==0 && flagQempty==1){
        
        //launch the A* kernel
        A_star<<<numThreads,numBlocks>>>(D_offset,D_edges,D_weights,D_hx,D_parent,D_Kqueue,D_Kq_size,D_lock,D_nextFlag,D_flagEnd,N,edge_size,K,endNode);

        //launch the kernel to setNV from flags
        cudaDeviceSynchronize();
        
        setNV<<<numThreads,new_Blocks>>>(D_nextFlag,D_nextV,D_nV_size,N);

        cudaDeviceSynchronize();
        
        //lauch kernel to insert in proirity q
        insertPQ<<<numThreads,numBlocks>>>(D_Kqueue,D_Kq_size,D_nextV,D_nV_size,K,N);

        //cpy flagend and flagEmpty
        cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(H_Kq_size,D_Kq_size, sizeof(int)*K,cudaMemcpyDeviceToHost);
        cudaMemcpy(D_nextFlag,H_nextFlag,sizeof(int)*N,cudaMemcpyHostToDevice);
        cudaMemcpy(D_nV_size,a0,sizeof(int),cudaMemcpyHostToDevice);

        cudaMemcpy(H_Kqueue,D_Kqueue, sizeof(int)*N,cudaMemcpyDeviceToHost);
       
        flagQempty = 0;
        for(int i=0;i<K;i++){
            // printf("size: %d\n",H_Kq_size[i]);
            if(H_Kq_size[i]>0)
                flagQempty=1;
        }
        
        // printf("fl; %d\n",*H_flagEnd);
        j++;
    }

//    cudaMemcpy(H_cx,D_cx, sizeof(int)*N,cudaMemcpyDeviceToHost);
    // cudaMemcpy(H_cx,Cx, sizeof(int)*N,cudaMemcpyDeviceToHost);
    cudaMemcpy(H_Kqueue,D_Kqueue, sizeof(int)*N,cudaMemcpyDeviceToHost);
    cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost);

    // printf("%d-%d cost is: %d\n",startNode,endNode,H_cx[endNode]);

    printf("path: %d ",endNode);
    int pr = H_parent[endNode];
    while(pr!=-1){
        printf("%d ",pr);
        pr = H_parent[pr];
    }
    printf("\n");
    
    return 0;

}
