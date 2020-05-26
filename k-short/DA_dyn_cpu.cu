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
#include <algorithm>


/***all macros**/
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

/***function declarations***/
void insertDiff(unordered_map< unsigned int, Node*>& Graph,int a,int b,unsigned int c);

void createDiffGraph(int N,unordered_map<unsigned int,Node*>& Graph,
    int* diffOff,int* diffEdges,unsigned int* diffWeight );

void removeDelEdges(int u,int v,int* offset,int* edges,int N,int E,int* rev_offset,int* rev_edges,int& del_size);

void mergeDiff(int* offset,int* edges,unsigned int* weight,int N,int E,
    int* diff_offset, int* diff_edges,unsigned int* diff_weight,int insert_size,int del_size,
    int* mOffset,int* mEdges,unsigned int* mWeight);

void check_del_path(int u, int v,vector<int> Path, bool& flag);
void check_cycle(int N,int* parent);

void computeTime(float& time,cudaEvent_t start, cudaEvent_t stop);

void check_remove_cycle(int* nodes,int* size,int* rev_offset,int* rev_edges,unsigned int* rev_weight,
                        int* parent,int* parent_old,int* Cx,int* Hx,int N,int E,
                        int* rev_diff_offset,int* rev_diff_edges,unsigned int* rev_diff_weight,int dE);

/**** device Code *******/

// __device__ volatile int Cx[MAX_NODE];
__device__ volatile int PQ[MAX_NODE];


//K in parallel
__global__ void extractMin(int* PQ_size, int* expandNodes,int* expandNodes_size,int* Cx,int* openList,int N,int K){
    
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
__global__ void A_star_expand(int* off,int* edge,unsigned int* W,int* Hx,int* parent,volatile int* Cx,
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagfound,int* openList,
    int N,int E, int K,int dest,int* nVFlag,int* PQ_size,
    int flagDiff,int* diff_off,int* diff_edge,unsigned int* diff_weight,int dE ){
       
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
                    if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
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
            
        }
        //end diff
        
    }//end 

}


//K in parallel -- O(N)
__global__ void keepHeapPQ(int* PQ_size,int* Cx,int N,int K){
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
        if(nextFlag[id]==1){
            int index = atomicAdd(nvSize,1);
            nextV[index]=id;
        }
    }
}


//for K in parallel
__global__ void insertPQ(int* PQS,int* nextV,int* nVsize,int* Cx,int K,int N,int* openList){
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
__global__ void checkMIN(int* PQ_size,int* flagEnd,int* Cx,int dest,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < K && PQ_size[id] > 0 ){
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];
        //check if atleast one min, dont end the a*
        if( Cx[node] < Cx[dest] ){
            atomicAnd(flagEnd,0);
        }
    }
}


__global__ void propogateDel(int* delEdgesV,int delEdge, volatile int* Cx,
                int* rev_offset,int* rev_edges,unsigned int* rev_weight,int N,int E,
                int* Hx,volatile int* parent,int* parent_old,int* addFlag,
                int* rev_diff_offset,int* rev_diff_edges,unsigned int* rev_diff_weight,int dE){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;

    if(id<delEdge){
        int node = delEdgesV[id];
        //check for the parent and add to nextflag and update the cost
        

        int start = rev_offset[node];
        int end = E;
        if(node!=N-1)
            end = rev_offset[node+1];
        
        //no parent
        // write in parent read always from old_parent
        parent[node] = -1;
        Cx[node]=INT_MAX;
        addFlag[node]=1;

        int cost = INT_MAX;
        int opt_parent = -1;

        //if any parent can change the cost 
        while(start< end){
            int p = rev_edges[start];
            
            //del edges
            if(p<0 || p==node){
                start++;
                continue;
            }
            
            int weight = rev_weight[start];
            int flag_cycle = false;
            
            //check parent doesn't contain node
            int ancestor = parent_old[p];
    
            while(ancestor>0){
                if(ancestor==node){
                    flag_cycle = true;
                    break;
                }
                ancestor = parent_old[ancestor];
                
            }
            
            
            //no need to lock only single parent so only one node in array so one node per thread
            if(!flag_cycle && Cx[p]!=INT_MAX && cost > (Cx[p]-Hx[p])+weight+Hx[node] ){
                cost = (Cx[p]-Hx[p] )+weight+Hx[node];
                opt_parent = p;
            }

            start++;
        }

        start = rev_diff_offset[node];
        end = dE;
        if(node!=N-1)
            end = rev_diff_offset[node+1];
        
        while(start< end){
            int p = rev_diff_edges[start];
            
            //del edges
            if(p<0 || p==node){
                start++;
                continue;
            }
            
            int weight = rev_diff_weight[start];
            int flag_cycle = false;
            
            //check parent doesn't contain node
            int ancestor = parent_old[p];
            
            while(ancestor!=-1){
                if(ancestor==node){
                    flag_cycle = true;
                    break;
                }
                ancestor = parent_old[ancestor];
                
            }

            //no need to lock only single parent so only one node in array so one node per thread
            if(!flag_cycle && Cx[p]!=INT_MAX && cost > (Cx[p]-Hx[p])+weight+Hx[node] ){
                cost = (Cx[p]-Hx[p] )+weight+Hx[node];
                opt_parent = p;
            }

            start++;
        }

        //write here
        if(cost!=INT_MAX){
            Cx[node]=cost;
            parent[node]=opt_parent;
        }

    }

}

//add inserted edges to propogate
__global__ void propogateAdd(int* diff_off, int* diff_edges,unsigned int* diff_W,int* Hx,int* addFlag,
            volatile int* Cx,int* lock, int* parent, int* parent_old, int N, int dE){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < N){
        int node = id;
        
        int start = diff_off[node];
        int end = dE;
        if(node!=N-1)
            end = diff_off[node+1];
        
        while(start < end ){
            int child = diff_edges[start];
            
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
                    bool flag_cycle = false;

                    int ancestor = node;
                    while(ancestor > 0){
                        if(ancestor==child){
                            flag_cycle = true;
                            break;
                        }
                        ancestor = parent_old[ancestor];
                        
                    }
                   
                    if(!flag_cycle && Cx[node] != INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
                        
                        Cx[child] =  (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child];
                
                        parent[child] = node;
                        __threadfence();

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


__global__ void insert_propagate(int* nodes, int* size, int* off, int* edge,unsigned int* W,int* Hx,
            int N,int E,volatile int* Cx,int* lock, int* parent,int* addFlag,
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
            bool leaveLoop = false;
                
            while(!leaveLoop){

                if(atomicExch(&lock[child],1)==0){
                    
                    if(Cx[node]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }

                    leaveLoop = true;    
                    atomicExch(&lock[child],0);
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
            bool leaveLoop = false;

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section

                    if(Cx[node]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
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


__global__ void delete_propagate(int* nodes, int* size, int* off, int* edge,unsigned int* W,int* Hx,
                    int N,int E,volatile int* Cx,int* lock, int* parent,int* parent_old,int* addFlag,
                    int* diff_off,int* diff_edge,unsigned int* diff_W,int dE,
                    int* rev_offset,int* rev_edges,unsigned int* rev_weight,
                    int* rev_diff_offset,int* rev_diff_edges,unsigned int* rev_diff_weight){

    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < *size){
        int node = nodes[id];
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start < end ){
            int child = edge[start];
            if(child<0){
                start++;
                continue;
            }

            bool leaveLoop = false;
                
            while(!leaveLoop){

                if(atomicExch(&lock[child],1)==0){
                    if(Cx[node]!=INT_MAX && Cx[child]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }
                    else 
                    if( (Cx[node]==INT_MAX && parent[child]==node ) || ( parent[child]==node && (Cx[child] < Cx[node] - Hx[node]+ W[start]+ Hx[child]) )  ){
                        //use back edges
                        int rstart = rev_offset[child];
                        int rend = E;
                        if(child!=N-1)
                            rend = rev_offset[child+1];
                        
                        //there is always one parent that is node.
                        Cx[child] = INT_MAX;
                        parent[child]=-1;

                        while(rstart < rend){
                            int p = rev_edges[rstart]; 
                            if(p<0 || p == child){
                                rstart++;
                                continue;
                            }

                            int weight = rev_weight[rstart];
                            bool flag_cycle = false;
                            
                            //check parent doesn't contain child
                        
                            int ancestor = parent_old[p];
                            
                            while(ancestor > 0){
                                if(ancestor==child){
                                    flag_cycle = true;
                                    break;
                                }
                                ancestor = parent_old[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }
                        
                        rstart =  rev_diff_offset[child];
                        rend = dE;
                        if(child!=N-1)
                            rend = rev_diff_offset[child+1];
                    
                        while(rstart < rend){
                            int p = rev_diff_edges[rstart]; 
                            
                            if(p<0 || p==child){
                                rstart++;
                                continue;
                            }

                            int weight = rev_diff_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent_old[p];
                            while(ancestor!=-1){
                                if(ancestor==child){
                                    flag_cycle = true;
                                    break;
                                }
                                    
                                ancestor = parent_old[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }

                        addFlag[child]=1;
                    }
                        

                    leaveLoop = true;
    
                    atomicExch(&lock[child],0);

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
                    if(Cx[node]!=INT_MAX && Cx[child]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }
                    else 
                    if((Cx[node]==INT_MAX && parent[child]==node )|| ( parent[child]==node && (Cx[child] < Cx[node] - Hx[node]+ diff_W[start]+ Hx[child]) )  ){
                        //use back edges
                        int rstart = rev_offset[child];
                        int rend = E;
                        if(child!=N-1)
                            rend = rev_offset[child+1];
                        
                        //there is always one parent that is node.
                        Cx[child] = INT_MAX;
                        parent[child]=-1;
 
                        while(rstart < rend){
                            int p = rev_edges[rstart]; 
                            
                            if(p<0 || p ==child){
                                rstart++;
                                continue;
                            }
 
                            int weight = rev_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent_old[p];
                            while(ancestor!=-1){
                                if(ancestor==child)
                                    flag_cycle = true;
                                ancestor = parent_old[ancestor];
                            }
                            
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }
 
                        rstart =  rev_diff_offset[child];
                        rend = dE;
                        if(child!=N-1)
                            rend = rev_diff_offset[child+1];
                    
                        while(rstart < rend){
                            int p = rev_diff_edges[rstart]; 
                            
                            if(p<0 || p==child){
                                rstart++;
                                continue;
                            }
 
                            int weight = rev_diff_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent_old[p];
                            while(ancestor!=-1){
                                if(ancestor==child){
                                 flag_cycle = true;
                                 break;
                                }
                                    
                                ancestor = parent_old[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }
 
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
__global__ void insertDest(int* PQ_size,int* Cx,int dest,int* openList){
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

__global__ void getCx(int* Cx,int dest,int* val){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id==0){
        *val = Cx[dest];
    }
}



/**** main function ****/
int main(){

    //the K PQ
    int K ;
    printf("Number of pq:");
    scanf("%d\n",&K);

    printf("source, destination:");
    int startNode,endNode;
    scanf("%d %d",&startNode,&endNode);

    int num_paths;
    printf("how many paths:");
    scanf("%d",&num_paths);

    FILE* fgraph = fopen("graph.txt","r");
    FILE* fgraph_rev = fopen("graph_op.txt","r");
    int N,E;
    fscanf(fgraph_rev,"%d %d\n",&N,&E);
    fscanf(fgraph,"%d %d\n",&N,&E);
    

    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*E);
    unsigned int* H_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);
    int* H_hx = (int*)malloc(sizeof(int)*N);
    int* H_cx = (int*)malloc(sizeof(int)*N);
    int* H_parent = (int*)malloc(sizeof(int)*N);
    int* H_parent_old = (int*)malloc(sizeof(int)*N);


    int* H_PQ = (int*)malloc(sizeof(int)*N);
    int* H_openList = (int*)malloc(sizeof(int)*N);
    int* H_PQ_size = (int*)malloc(sizeof(int)*K);



    //for reverse graph
    int* H_rev_edges = (int*)malloc(sizeof(int)*E);
    int* H_rev_offset = (int*)malloc(sizeof(int)*N);
    unsigned int* H_rev_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);


    //for cost of endNode
    int* H_dest_cost = (int*)malloc(sizeof(int));
    
    memset(H_PQ_size,0,sizeof(int)*K);
    memset(H_openList,-1,sizeof(int)*N);

    //init cx
    for(int i=0;i<N;i++){
        H_cx[i]=INT_MAX;
        H_parent[i]=-1;
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%d",&H_edges[i]);
        fscanf(fgraph_rev,"%d",&H_rev_edges[i]);
    }

    for(int i=0;i<N;i++){
        fscanf(fgraph,"%d",&H_offset[i]);
        fscanf(fgraph_rev,"%d",&H_rev_offset[i]);
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%u",&H_weight[i]);
        fscanf(fgraph_rev,"%u",&H_rev_weight[i]);
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
    fclose(fgraph_rev);
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

    //insert startNode in PQ[0]
    H_cx[startNode]=H_hx[startNode];
    H_PQ[0]=startNode;
    H_PQ_size[0]=1;
    H_openList[startNode]=0;


    //create events to record runtime
    float run_time = 0;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //graph struture
    int* D_offset;
    int* D_edges ;
    unsigned int* D_weight;
    int* D_hx;
    int* D_parent;

    //for reading the ancessostor to avoid lock for write after read.
    int* D_parent_old;
    
    //Priority queue size
    int* D_PQ_size;

    //CX
    int* D_Cx;

    //flag if in openList(contains which PQ)
    int* D_openList;
    //lock for nodes
    int* D_lock;

    //Diff structure
    int* D_diff_edges;
    int* D_diff_offset;
    unsigned int* D_diff_weight;

    //reverse graph
    int* D_rev_edges;
    int* D_rev_offset;
    unsigned int* D_rev_weight;

    //reverse diff
    int* D_rev_diff_offset;
    int* D_rev_diff_edges;
    unsigned int* D_rev_diff_weight;

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
    
    //list of nodes v of deleted edges u->
    int* D_delEdgesV;

    gpuErrchk ( cudaMalloc(&D_offset,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_weight,sizeof(unsigned int)*E) );
    gpuErrchk ( cudaMalloc(&D_hx,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_parent,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_parent_old,sizeof(int)*N) );
   
    gpuErrchk ( cudaMalloc(&D_Cx,sizeof(int)*N) );

    gpuErrchk ( cudaMalloc(&D_PQ_size,sizeof(int)*K) );
    
    gpuErrchk ( cudaMalloc(&D_openList,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_lock,sizeof(int)*N) );

    gpuErrchk ( cudaMalloc(&D_dest_cost,sizeof(int)) );

    //rev graph
    gpuErrchk ( cudaMalloc(&D_rev_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_rev_offset,sizeof(int)*N ) );   
    gpuErrchk ( cudaMalloc(&D_rev_weight,sizeof(unsigned int)*E) );


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

    gpuErrchk ( cudaMemcpy(D_Cx,H_cx,sizeof(int)*N,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpyToSymbol(PQ,H_PQ, sizeof(int)*N, 0, cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_flagfound,H_flagfound,sizeof(int),cudaMemcpyHostToDevice) );
   
    gpuErrchk ( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_expandNodes_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );

    //reverse graph
    gpuErrchk ( cudaMemcpy(D_rev_offset,H_rev_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_rev_edges,H_rev_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(D_rev_weight,H_rev_weight,sizeof(unsigned int)*E,cudaMemcpyHostToDevice) );


    
    gpuErrchk ( cudaMemset(D_lock,0,sizeof(int)*N) );


    int flag_PQ_not_empty = 0;
    for(int i=0;i<K;i++){
        if(H_PQ_size[i]>0)
            flag_PQ_not_empty=1;
    }


    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int N_numBlocks = (N+numThreads-1)/numThreads;

    if(DEBUG)
        printf("[INFO] A* started\n");


    cudaEventRecord(start);
   
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

    getCx<<<1,1>>>(D_Cx,endNode,D_dest_cost);
    gpuErrchk( cudaMemcpy(H_dest_cost,D_dest_cost, sizeof(int),cudaMemcpyDeviceToHost) );

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    computeTime(run_time,start,stop);

 
    gpuErrchk( cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );


    vector<int> Path;
    printf("[OUT] Cost: %d\n",*H_dest_cost);
    printf("[OUT] Path(in reverse): ");
    if(*H_dest_cost!=INT_MAX){
        int p = endNode;
        while(H_parent[p]!=-1){
            printf("%d ",p); 
            Path.push_back(p); 
            p = H_parent[p];
        }
        printf("%d\n",p);
        Path.push_back(p);
    }
    else{
        printf("not found\n");
    }

    reverse(Path.begin(),Path.end());



    ///////////////////////////////////////////////
    // A star complete //

    int line;
    int update_count = 0;
    for(int p=0;p<num_paths;p++){
        
        //list of nodes v of deleted edges u->v
        int* H_delEdgesV = (int*)malloc(sizeof(int)*E);
        gpuErrchk ( cudaMalloc(&D_delEdgesV,sizeof(int)*E) );

        unordered_map<unsigned int,Node*> Graph;
        unordered_map<unsigned int,Node*> rev_Graph;

        vector<pair<int,int>>deleted_edges;

        bool flag_do_a_star = false;
        
        int insertEdge=0, delEdge=0;
        int delEdgesV_size = 0;                    //v whose cost can change due to deletion

        // for(int i=0;i<line;i++){
        //     int flag;
        //     int u,v;
        //     unsigned int w;
        //     fscanf(fdiff,"%d %d %d %u\n",&flag,&u,&v,&w);
        //     if(flag==1){
        //         insertDiff(Graph,u,v,w);
        //         insertDiff(rev_Graph,v,u,w);
        //         insertEdge++;
        //     }
        //     else if(flag==0){
                
        //         deleted_edges.push_back(pair<int,int>(u,v));               
        //     }
            
        // }

        /* remove last edge from graph */
        if(Path.size()==0){
            break;
        }

        deleted_edges.push_back(pair<int,int>(Path[Path.size()-2],Path[Path.size()-1])); 


        // insertEdge is insertion size
        //for diff
        int* H_diff_edges = (int*)malloc(sizeof(int)*insertEdge);
        int* H_diff_offset = (int*)malloc(sizeof(int)*N);
        unsigned int* H_diff_weight = (unsigned int*)malloc(sizeof(unsigned int)*insertEdge);

        //diff for revrse graph
        int* H_rev_diff_edges = (int*)malloc(sizeof(int)*insertEdge);
        int* H_rev_diff_offset = (int*)malloc(sizeof(int)*N);
        unsigned int* H_rev_diff_weight = (unsigned int*)malloc(sizeof(unsigned int)*insertEdge);

        //diff csr
        gpuErrchk ( cudaMalloc(&D_diff_edges,sizeof(int)*insertEdge) );
        gpuErrchk ( cudaMalloc(&D_diff_offset,sizeof(int)*(N+1) ) );   //coz
        gpuErrchk ( cudaMalloc(&D_diff_weight,sizeof(unsigned int)*insertEdge) );
        

        //rev diff graph
        gpuErrchk ( cudaMalloc(&D_rev_diff_edges,sizeof(int)*insertEdge) );
        gpuErrchk ( cudaMalloc(&D_rev_diff_offset,sizeof(int)*(N+1) ) );   
        gpuErrchk ( cudaMalloc(&D_rev_diff_weight,sizeof(unsigned int)*insertEdge) );


        //reset offset to 0 ..ie no nodes
        memset(H_diff_offset,0,sizeof(int)*N);
        memset(H_rev_diff_offset,0,sizeof(int)*N);

        if(1)
            printf("[INFO](%d) insertion:%d\n",update_count,insertEdge);

        createDiffGraph(N,Graph,H_diff_offset,H_diff_edges,H_diff_weight);
        createDiffGraph(N,rev_Graph,H_rev_diff_offset,H_rev_diff_edges,H_rev_diff_weight);
        
        //TODO free the graphs

        //diff graph
        gpuErrchk ( cudaMemcpy(D_diff_edges,H_diff_edges,sizeof(int)*insertEdge,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_diff_offset,H_diff_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_diff_weight,H_diff_weight,sizeof(unsigned int)*insertEdge,cudaMemcpyHostToDevice) );

        //rev diff graph
        gpuErrchk ( cudaMemcpy(D_rev_diff_edges,H_rev_diff_edges,sizeof(int)*insertEdge,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_diff_offset,H_rev_diff_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_diff_weight,H_rev_diff_weight,sizeof(unsigned int)*insertEdge,cudaMemcpyHostToDevice) );

        //reset D_nV flag
        gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );


        //do insertion first
        if(DEBUG)
            printf("[INFO] starting computing cost for inserions\n");

    
        gpuErrchk( cudaMemcpy(D_parent_old,D_parent,sizeof(int)*N,cudaMemcpyDeviceToDevice) );
        
        cudaEventRecord(start);
        //N parallel
        propogateAdd<<<N_numBlocks,numThreads>>>(D_diff_offset, D_diff_edges,D_diff_weight,D_hx,D_nVFlag,
            D_Cx,D_lock,D_parent,D_parent_old,N,insertEdge);
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        computeTime(run_time,start,stop);

        
        gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );

        
        //gen from flag D_nV
        cudaEventRecord(start);
        setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        computeTime(run_time,start,stop);

        //copy back
        gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );

        //reset nV flags
        gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );
 
        if(DEBUG)
             printf("[INFO] starting propogation for insertions\n");
 
        while(*H_nV_size > 0){ 

            numBlocks = (*H_nV_size+numThreads-1)/numThreads;
            cudaEventRecord(start);
            insert_propagate<<<numBlocks,numThreads>>>(D_nV,D_nV_size,D_offset,D_edges,D_weight,D_hx,
                N,E,D_Cx,D_lock,D_parent,D_nVFlag,
                D_diff_offset,D_diff_edges,D_diff_weight,insertEdge);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            computeTime(run_time,start,stop);
            
            //reset size=0
            gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );

            //gen from flag D_nV
            cudaEventRecord(start);
            setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            computeTime(run_time,start,stop);
            

            //copy back
            gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );

            //reset nV flags
            gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

        }

        getCx<<<1,1>>>(D_Cx,endNode,D_dest_cost);
        gpuErrchk( cudaMemcpy(H_dest_cost,D_dest_cost, sizeof(int),cudaMemcpyDeviceToHost) );
        
        //copy parent
        gpuErrchk( cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );

        Path.clear();
        if(*H_dest_cost!=INT_MAX){
            int p = endNode;
            while(H_parent[p]!=-1){
                Path.push_back(p);
                p = H_parent[p];
            }
            Path.push_back(p);
        }
        //reverse the path to get from source to end
        reverse(Path.begin(),Path.end());

        //start computation for deletion
        for(int j=0;j<deleted_edges.size();j++){
            int u,v;
            u= deleted_edges[j].first;
            v= deleted_edges[j].second;

            //check id del edges in optimal path.
            check_del_path(u,v,Path,flag_do_a_star);

            int old = delEdge;
            //if deleted adds to delEdge
            removeDelEdges(u,v,H_offset,H_edges,N,E,H_rev_offset,H_rev_edges,delEdge);
            //add to list only if its cost changes due to this deletion
            if(H_parent[v]==u && delEdge > old ){
                H_delEdgesV[delEdgesV_size]=v;
                delEdgesV_size++;
            }

        }

        //free memory
        deleted_edges.clear();
        deleted_edges.shrink_to_fit();      //may not but gives a hint


        //deleted edges
        gpuErrchk ( cudaMemcpy(D_edges,H_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_edges,H_rev_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_delEdgesV,H_delEdgesV,sizeof(int)*E,cudaMemcpyHostToDevice) );


        //remove check
        // check_cycle(N,H_parent);

        if(1)
            printf("[INFO](%d) deletion:%d, eff del:%d\n",update_count,delEdge,delEdgesV_size);


        //add del
        if(delEdgesV_size>0){
            if(DEBUG)
                printf("[INFO] Starting computing cost for deletions\n");

            //old parent to check cycle
            gpuErrchk( cudaMemcpy(D_parent_old,D_parent,sizeof(int)*N,cudaMemcpyDeviceToDevice) );

            int numBlocks_del = ( delEdgesV_size + numThreads -1)/numThreads;
            
            cudaEventRecord(start);
            propogateDel<<<numBlocks_del,numThreads>>>(D_delEdgesV,delEdgesV_size,D_Cx,
                D_rev_offset,D_rev_edges,D_rev_weight,N,E,
                D_hx,D_parent,D_parent_old,D_nVFlag,
                D_rev_diff_offset,D_rev_diff_edges,D_rev_diff_weight,insertEdge);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();

            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            computeTime(run_time,start,stop);

            //check
            gpuErrchk( cudaMemcpy(H_parent,D_parent,sizeof(int)*N,cudaMemcpyDeviceToHost) );
            gpuErrchk( cudaMemcpy(H_parent_old,D_parent_old,sizeof(int)*N,cudaMemcpyDeviceToHost) );
            gpuErrchk( cudaMemcpy(H_cx,D_Cx,sizeof(int)*N,cudaMemcpyDeviceToHost) );

            //remove cycles
            check_remove_cycle(H_delEdgesV,&delEdgesV_size,H_rev_offset,H_rev_edges,H_rev_weight,
                H_parent,H_parent_old,H_cx,H_hx,N,E,
                H_rev_diff_offset,H_rev_diff_edges,H_rev_diff_weight,insertEdge);

            gpuErrchk( cudaMemcpy(D_Cx,H_cx,sizeof(int)*N,cudaMemcpyHostToDevice) );
            gpuErrchk( cudaMemcpy(D_parent,H_parent,sizeof(int)*N,cudaMemcpyHostToDevice) );

            // check_cycle(N,H_parent);

            if(DEBUG)
                printf("[INFO] starting propogation for deletions\n");

            //make size =0
            gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
        
            //gen from flag D_nV
            cudaEventRecord(start);
            setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            computeTime(run_time,start,stop);

            //copy back
            gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );

            //reset nV flags
            gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );
            

            while(*H_nV_size > 0){ 

                numBlocks = (*H_nV_size+numThreads-1)/numThreads;

               //printf("delsize: %d\n",*H_nV_size);

                //old parent to check cycle and remove locking on parent
                gpuErrchk( cudaMemcpy(D_parent_old,D_parent,sizeof(int)*N,cudaMemcpyDeviceToDevice) );
               
                cudaEventRecord(start);
                delete_propagate<<<numBlocks,numThreads>>>(D_nV,D_nV_size,D_offset,D_edges,D_weight,D_hx,
                    N,E,D_Cx,D_lock,D_parent,D_parent_old,D_nVFlag,
                    D_diff_offset,D_diff_edges,D_diff_weight,insertEdge,
                    D_rev_offset,D_rev_edges,D_rev_weight,
                    D_rev_diff_offset,D_rev_diff_edges,D_rev_diff_weight);
                
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);


                //reset size=0
                gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
    
                //gen from flag D_nV
                cudaEventRecord(start);
                setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
                
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);
                
    
                //copy back
                gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );
    
                //reset nV flags
                gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

                //check cycle
                gpuErrchk( cudaMemcpy(H_parent,D_parent,sizeof(int)*N,cudaMemcpyDeviceToHost) );
                gpuErrchk( cudaMemcpy(H_parent_old,D_parent_old,sizeof(int)*N,cudaMemcpyDeviceToHost) );
                gpuErrchk( cudaMemcpy(H_cx,D_Cx,sizeof(int)*N,cudaMemcpyDeviceToHost) );
                gpuErrchk( cudaMemcpy(H_nV,D_nV,sizeof(int)*N,cudaMemcpyDeviceToHost) );


                //remove cycles
                check_remove_cycle(H_nV,H_nV_size,H_rev_offset,H_rev_edges,H_rev_weight,
                    H_parent,H_parent_old,H_cx,H_hx,N,E,
                    H_rev_diff_offset,H_rev_diff_edges,H_rev_diff_weight,insertEdge);

                gpuErrchk( cudaMemcpy(D_Cx,H_cx,sizeof(int)*N,cudaMemcpyHostToDevice) );
                gpuErrchk( cudaMemcpy(D_parent,H_parent,sizeof(int)*N,cudaMemcpyHostToDevice) );
    
            }
        
        }

        if(DEBUG)
            printf("[INFO] updating priority queue\n");

        //propogate complete do normal A*
        numBlocks = (K+numThreads-1)/numThreads;

        //update PQ after propogate
        cudaEventRecord(start);
        keepHeapPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_Cx,N,K);
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        computeTime(run_time,start,stop);

        //check if there is node cost in PQ less than dest
        *H_flagEnd = 1;
        gpuErrchk( cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice) );
    
        cudaEventRecord(start);
        checkMIN<<< numBlocks,numThreads >>>(D_PQ_size,D_flagEnd,D_Cx,endNode,N,K);
        
        gpuErrchk( cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        computeTime(run_time,start,stop);

        gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );
        

        //here flag end represents from above that there is a node with cost lesser     
        if(*H_flagEnd==0 && flag_do_a_star){

            printf("[INFO] doing a* after propogation\n");

            cudaEventRecord(start);
            insertDest<<<1,1>>>(D_PQ_size,D_Cx,endNode,D_openList);
            gpuErrchk(cudaPeekAtLastError() );
            
            cudaDeviceSynchronize();
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            computeTime(run_time,start,stop);


            flag_PQ_not_empty = 0;
            for(int i=0;i<K;i++){
                if(H_PQ_size[i]>0)
                    flag_PQ_not_empty=1;
            }
    
            //reset flags
            *H_flagEnd = 0;
            *H_flagfound = 0;
            gpuErrchk ( cudaMemcpy(D_flagfound,H_flagfound,sizeof(int),cudaMemcpyHostToDevice) );
            gpuErrchk( cudaMemcpy(D_nV_size,H_a0,sizeof(int),cudaMemcpyHostToDevice) );
            gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );
    
            //DO A* initailly on whole graph
            while(*H_flagEnd==0 && flag_PQ_not_empty==1){
                
                //extract min
                cudaEventRecord(start);
                extractMin<<<numBlocks,numThreads>>>(D_PQ_size,D_expandNodes,D_expandNodes_size,D_Cx,D_openList,N,K);
                
                gpuErrchk(cudaPeekAtLastError() );
    
                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);

    
                cudaEventRecord(start);
                A_star_expand<<<numBlocks,numThreads>>>(D_offset,D_edges,D_weight,D_hx,D_parent,D_Cx,
                    D_expandNodes,D_expandNodes_size, D_lock ,D_flagfound,D_openList,
                    N,E,K,endNode,D_nVFlag,D_PQ_size,
                    true,D_diff_offset,D_diff_edges,D_diff_weight,insertEdge);
                
                gpuErrchk(cudaPeekAtLastError() );
                
                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);


                cudaEventRecord(start);
                keepHeapPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_Cx,N,K);
                gpuErrchk(cudaPeekAtLastError() );

                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);

                
                //gen from flag D_nV
                //for N in parallel
                cudaEventRecord(start);
                setNV<<<N_numBlocks,numThreads>>>(D_nVFlag,D_nV,D_nV_size,N);
                
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);

                
                cudaEventRecord(start);
                insertPQ<<<numBlocks,numThreads>>>(D_PQ_size,D_nV,D_nV_size,D_Cx,K,N,D_openList);
                
                gpuErrchk(cudaPeekAtLastError() );
                cudaDeviceSynchronize();
                cudaEventRecord(stop);
                cudaEventSynchronize(stop);
                computeTime(run_time,start,stop);

            
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

                    cudaEventRecord(start);
                    checkMIN<<< numBlocks,numThreads >>>(D_PQ_size,D_flagEnd,D_Cx,endNode,N,K);
                    
                    gpuErrchk( cudaPeekAtLastError() );
                    cudaDeviceSynchronize();
                    cudaEventRecord(stop);
                    cudaEventSynchronize(stop);
                    computeTime(run_time,start,stop);
    

                    gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );
                
                }
        
            }
        }
    
        cudaEventRecord(start);
        getCx<<<1,1>>>(D_Cx,endNode,D_dest_cost);
        
        cudaDeviceSynchronize();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        computeTime(run_time,start,stop);
    
        gpuErrchk( cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );


        // found or not found based on Cx
        gpuErrchk( cudaMemcpy(H_dest_cost,D_dest_cost, sizeof(int),cudaMemcpyDeviceToHost) );

        //remove old path
        Path.clear();
        printf("[OUT] Cost: %d\n",*H_dest_cost);
        printf("[OUT] Path(in reverse): ");
        if(*H_dest_cost!=INT_MAX){
            int p = endNode;
            while(H_parent[p]!=-1){
                printf("%d ",p);
                Path.push_back(p);
                p = H_parent[p];
            }
            Path.push_back(p);
            printf("%d\n",p);
        }
        else{
            printf("not found\n");
        }

        //reverse the path to get from source to end
        reverse(Path.begin(),Path.end());

        
        //merge graph
        int* H_offset_new,*H_edges_new;
        unsigned int* H_weight_new;
        
        int E_new = E + insertEdge - delEdge;

        H_offset_new = (int*)malloc(sizeof(int)*N);
        H_edges_new =  (int*)malloc(sizeof(int)*E_new);
        H_weight_new = (unsigned int*)malloc(sizeof(unsigned int)*E_new);

        mergeDiff(H_offset,H_edges,H_weight,N,E,
            H_diff_offset,H_diff_edges,H_diff_weight,insertEdge,delEdge,
            H_offset_new,H_edges_new,H_weight_new);

        //free pointer
        free(H_offset);
        free(H_edges);
        free(H_weight);
        free(H_diff_offset);
        free(H_diff_edges);
        free(H_diff_weight);

        H_offset = H_offset_new;
        H_edges = H_edges_new;
        H_weight = H_weight_new;
        
        //cudaFree and cpy
        cudaFree(D_edges);
        cudaFree(D_weight);
        cudaFree(D_diff_edges);
        cudaFree(D_diff_offset);
        cudaFree(D_diff_weight);

        gpuErrchk ( cudaMalloc(&D_edges,sizeof(int)*E_new) );
        gpuErrchk ( cudaMalloc(&D_weight,sizeof(unsigned int)*E_new) );
   
        gpuErrchk ( cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_edges,H_edges,sizeof(int)*E_new,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_weight,H_weight,sizeof(unsigned int)*E_new,cudaMemcpyHostToDevice) );

        //merge rev graph
        int* H_rev_offset_new,*H_rev_edges_new;
        unsigned int* H_rev_weight_new;

        H_rev_offset_new = (int*)malloc(sizeof(int)*N);
        H_rev_edges_new =  (int*)malloc(sizeof(int)*E_new);
        H_rev_weight_new = (unsigned int*)malloc(sizeof(unsigned int)*E_new);
        
        mergeDiff(H_rev_offset,H_rev_edges,H_rev_weight,N,E,
            H_rev_diff_offset,H_rev_diff_edges,H_rev_diff_weight,insertEdge,delEdge,
            H_rev_offset_new,H_rev_edges_new,H_rev_weight_new);
        
        free(H_rev_offset);
        free(H_rev_edges);
        free(H_rev_weight);
        free(H_rev_diff_offset);
        free(H_rev_diff_edges);
        free(H_rev_diff_weight);
        
        H_rev_offset = H_rev_offset_new;
        H_rev_edges = H_rev_edges_new;
        H_rev_weight = H_rev_weight_new;

        //cuda free and cpy
        cudaFree(D_rev_edges);
        cudaFree(D_rev_weight);
        cudaFree(D_rev_diff_edges);
        cudaFree(D_rev_diff_offset);
        cudaFree(D_rev_diff_weight);

        gpuErrchk ( cudaMalloc(&D_rev_edges,sizeof(int)*E_new) );
        gpuErrchk ( cudaMalloc(&D_rev_weight,sizeof(unsigned int)*E_new) );

        gpuErrchk ( cudaMemcpy(D_rev_offset,H_rev_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_edges,H_rev_edges,sizeof(int)*E_new,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_weight,H_rev_weight,sizeof(unsigned int)*E_new,cudaMemcpyHostToDevice) );

        //change E
        E = E_new;
        
        cudaFree(D_delEdgesV);
        free(H_delEdgesV);

        //inc
        update_count++;
        
    }

    printf("[INFO] update count: %d\n",update_count);
    printf("[INFO] RUNTIME: %f\n",run_time);

    //cuda free
    // free everything

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

            if(offindex < N ){
                diffOff[offindex] = k;
                offindex++;
            }
            
        }
        else{
            if(offindex < N ){
                diffOff[offindex] = k;
                offindex++;
            }
           
        }
    }

}


void removeDelEdges(int u,int v,int* offset,int* edges,int N,int E,int* rev_offset,int* rev_edges,int& del_size){
    int start = offset[u];
    int end =  E;
    bool flag_done = false;
    bool flag_done_rev = false;
    if(u!=N-1)
        end = offset[u+1];
    while(start<end){
        if( v == edges[start]){
            edges[start]=-1;
            flag_done = true;
            break;
        }
        start++;
    }

    start = rev_offset[v];
    end = E;
    if(v!=N-1)
        end = rev_offset[v+1];
    while(start < end){
        if(u == rev_edges[start]){
            rev_edges[start] = -1;
            flag_done_rev = true;
            break;
        }
        start++;
    }

    if(flag_done && flag_done_rev)
        del_size++;
    
    if( (flag_done && !flag_done_rev)|| (!flag_done && flag_done_rev) )
        printf("[ERROR] edge present in front ot back graph\n");

}


void check_del_path(int u, int v,vector<int> Path, bool& flag){
    vector<int> :: iterator itr;
    itr = find(Path.begin(),Path.end(),u);
    if(itr!=Path.end()){
        itr+=1;
        if(*itr == v)
            flag = true;
    }

}

void check_cycle(int N,int* parent){
    int flag = 0;
    for(int i=0;i<N;i++){
        
        vector<int> visited(N,0);
        int ancestor = parent[i];
        while(ancestor > 0){
            if(visited[ancestor]==1){
                printf("cycle at: %d, %d\n",i,ancestor); 
                flag =1;
                break;
            }

            visited[ancestor]=1;
            ancestor = parent[ancestor];
        }
        
    }
    if(flag==0)
        printf("no cycle\n");
}


void mergeDiff(int* offset,int* edges,unsigned int* weight,int N,int E,
            int* diff_offset, int* diff_edges,unsigned int* diff_weight,int insert_size,int del_size,
            int* mOffset,int* mEdges,unsigned int* mWeight){

  
    int E_new = E + insert_size - del_size;
    mOffset[0] = 0;
    int edegOffset= 0;
    for(int i=0;i<N;i++){

        int start = offset[i];
        int end = E;
        if(i!=N-1)
            end = offset[i+1];
        //int count = 0;
        while(start<end){
            int child = edges[start];
            if(child!=-1){
                mEdges[edegOffset] = child;
                mWeight[edegOffset] = weight[start];
                edegOffset++;
            }
            start++;
        }

        start = diff_offset[i];
        end =  insert_size;
        if(i!=N-1)
            end = diff_offset[i+1];
        while(start<end){
            int child = diff_edges[start];
            if(child!=-1){
                mEdges[edegOffset] = child;
                mWeight[edegOffset]= diff_weight[start];
                edegOffset++;
            }
            start++;
        }

        if(edegOffset > E_new){
            printf("ERROR: size %d::%d\n",E_new,edegOffset);
            break;
        }

        if(i!=N-1)
            mOffset[i+1]=edegOffset;
    }

}

void computeTime(float& time,cudaEvent_t start, cudaEvent_t stop){
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    time+= milliseconds;
    //printf("time:%f\n",milliseconds);
}


void check_remove_cycle(int* nodes,int* size,int* rev_offset,int* rev_edges,unsigned int* rev_weight,
            int* parent,int* parent_old,int* Cx,int* Hx, int N,int E,
            int* rev_diff_offset,int* rev_diff_edges,unsigned int* rev_diff_weight,int dE){

    for(int i=0;i<*size;i++){
        int node = nodes[i];
        bool cycle = false;

        vector<bool>visited(N,false);
        int ancestor = parent[node];
        while(ancestor > 0){
            if(ancestor==node){
                cycle  = true;
                break;
            }
            if(visited[ancestor]){
                break;
            }
            visited[ancestor]=true;
            ancestor = parent[ancestor];
        }

        if(cycle){
          // printf("cycle at %d->%d\n",node,parent[node]);
            int p_cycle = parent[node];

            int start = rev_offset[node];
            int end = E;
            if(node!=N-1)
                end = rev_offset[node+1];
            
            //no parent
            // write in parent read always from old_parent

            int cost = INT_MAX;
            int opt_parent = -1;

            //if any parent can change the cost 
            while(start< end){
                int p = rev_edges[start];
                
                //del edges
                if(p<0 || p == p_cycle || p==node){
                    start++;
                    continue;
                }
                
                int weight = rev_weight[start];
                int flag_cycle = false;
                
                //check parent doesn't contain node
                int ancestor = parent_old[p];
                while(ancestor>0){
                    if(ancestor==node){
                        flag_cycle = true;
                        break;
                    }
                    ancestor = parent_old[ancestor];
                }

                //no need to lock only single parent so only one node in array so one node per thread
                if(!flag_cycle && Cx[p]!=INT_MAX && cost > (Cx[p]-Hx[p])+weight+Hx[node] ){
                    cost = (Cx[p]-Hx[p] )+weight+Hx[node];
                    opt_parent = p;
                }

                start++;
            }

            start = rev_diff_offset[node];
            end = dE;
            if(node!=N-1)
                end = rev_diff_offset[node+1];
            
            while(start< end){
                int p = rev_diff_edges[start];
                
                //del edges
                if(p<0 || p== p_cycle || p== node){
                    start++;
                    continue;
                }
                
                int weight = rev_diff_weight[start];
                int flag_cycle = false;
                
                //check parent doesn't contain node
                int ancestor = parent_old[p];
               
                while(ancestor!=-1){
                    if(ancestor==node){
                        flag_cycle = true;
                        break;
                    }
                    ancestor = parent_old[ancestor];
                   
                }
               
                //no need to lock only single parent so only one node in array so one node per thread
                if(!flag_cycle && Cx[p]!=INT_MAX && cost > (Cx[p]-Hx[p])+weight+Hx[node] ){
                    cost = (Cx[p]-Hx[p] )+weight+Hx[node];
                    opt_parent = p;
                }

                start++;
            }

        
            Cx[node]=cost;
            parent[node]=opt_parent;
          // printf("opt p : %d\n",opt_parent);

        }
    
    }

}