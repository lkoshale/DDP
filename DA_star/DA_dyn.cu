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

void removeDelEdges(int u,int v,int* offset,int* edges,int N,int E,int* rev_offset,int* rev_edges);

void check_del_path(int u, int v,vector<int> Path, bool& flag);

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
        if(nextFlag[id]==1){
            int index = atomicAdd(nvSize,1);
            nextV[index]=id;
        }
    }
}


//for K in parallel
__global__ void insertPQ(int* PQS,int* nextV,int* nVsize,int K,int N,int* openList){
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
__global__ void checkMIN(int* PQ_size,int* flagEnd,int dest,int N,int K){
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


__global__ void propogateDel(int* delEdgesV,int delEdge,int* rev_offset,int* rev_edges,unsigned int* rev_weight,int N,int E,
                int* Hx,int* parent,int* lock,int* addFlag){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id<delEdge){
        int node = delEdgesV[id];
        //check for the parent and add to nextflag and update the cost
        int start = rev_offset[node];
        int end = E;
        if(node!=N-1)
            end = rev_offset[node+1];
        
        //no parent
        parent[node]=-1;
        Cx[node]=INT_MAX;
        addFlag[node]=1;

        //if any parent can change the cost 
        while(start< end){
            int p = rev_edges[start];
            
            //del edges
            if(p<0){
                start++;
                continue;
            }
            
            int weight = rev_weight[start];
            int flag_cycle = false;
            
            
            //check parent doesn't contain node
            int ancestor = parent[p];

            while(ancestor!=-1){
                if(ancestor==node)
                    flag_cycle = true;
                ancestor = parent[ancestor];
            }

            //no need to lock only single parent so only one node in array so one node per thread
            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[node] > (Cx[p]-Hx[p])+weight+Hx[node] ){
                Cx[node] = (Cx[p]-Hx[p] )+weight+Hx[node];
                parent[node] = p;
            }

            start++;
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

                    } else 
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
                            
                            if(p<0){
                                rstart++;
                                continue;
                            }

                            int weight = rev_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent[p];
                            while(ancestor!=-1){
                                if(ancestor==child)
                                    flag_cycle = true;
                                ancestor = parent[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }

                        //newly added backedges
                        rstart =  rev_diff_offset[child];
                        rend = dE;
                        if(child!=N-1)
                            rend = rev_diff_offset[child+1];
                    
                        while(rstart < rend){
                            int p = rev_diff_edges[rstart]; 
                            
                            if(p<0){
                                rstart++;
                                continue;
                            }

                            int weight = rev_diff_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent[p];
                            while(ancestor!=-1){
                                if(ancestor==child)
                                    flag_cycle = true;
                                ancestor = parent[ancestor];
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

                    if( Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }else 
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
                           
                           if(p<0){
                               rstart++;
                               continue;
                           }

                           int weight = rev_weight[rstart];
                           int flag_cycle = false;
                           
                           //check parent doesn't contain child
                           int ancestor = parent[p];
                           while(ancestor!=-1){
                               if(ancestor==child)
                                   flag_cycle = true;
                               ancestor = parent[ancestor];
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
                           
                           if(p<0){
                               rstart++;
                               continue;
                           }

                           int weight = rev_diff_weight[rstart];
                           int flag_cycle = false;
                           
                           //check parent doesn't contain child
                           int ancestor = parent[p];
                           while(ancestor!=-1){
                               if(ancestor==child)
                                   flag_cycle = true;
                               ancestor = parent[ancestor];
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


    int* H_PQ = (int*)malloc(sizeof(int)*N);
    int* H_openList = (int*)malloc(sizeof(int)*N);
    int* H_PQ_size = (int*)malloc(sizeof(int)*K);


    //for diff
    int* H_diff_edges = (int*)malloc(sizeof(int)*E);
    int* H_diff_offset = (int*)malloc(sizeof(int)*N);
    unsigned int* H_diff_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);

    //for reverse graph
    int* H_rev_edges = (int*)malloc(sizeof(int)*E);
    int* H_rev_offset = (int*)malloc(sizeof(int)*N);
    unsigned int* H_rev_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);
   
    //diff for revrse graph
    int* H_rev_diff_edges = (int*)malloc(sizeof(int)*E);
    int* H_rev_diff_offset = (int*)malloc(sizeof(int)*N);
    unsigned int* H_rev_diff_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);


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

    //required coz if many tries to add same in diff threads high low lower
    int* H_nVFlag = (int*)malloc(sizeof(int)*N);
    memset(H_nVFlag,-1,sizeof(int)*N);    

    //list of nodes v of deleted edges u->v
    int* H_delEdgesV = (int*)malloc(sizeof(int)*E);

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
   
    gpuErrchk ( cudaMalloc(&D_PQ_size,sizeof(int)*K) );
    
    gpuErrchk ( cudaMalloc(&D_openList,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&D_lock,sizeof(int)*N) );

    gpuErrchk ( cudaMalloc(&D_dest_cost,sizeof(int)) );
    gpuErrchk ( cudaMalloc(&D_delEdgesV,sizeof(int)*E) );


    //diff csr
    gpuErrchk ( cudaMalloc(&D_diff_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_diff_offset,sizeof(int)*(N+1) ) );   //coz
    gpuErrchk ( cudaMalloc(&D_diff_weight,sizeof(unsigned int)*E) );

    //rev graph
    gpuErrchk ( cudaMalloc(&D_rev_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_rev_offset,sizeof(int)*N ) );   
    gpuErrchk ( cudaMalloc(&D_rev_weight,sizeof(unsigned int)*E) );

    //rev diff graph
    gpuErrchk ( cudaMalloc(&D_rev_diff_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&D_rev_diff_offset,sizeof(int)*N ) );   
    gpuErrchk ( cudaMalloc(&D_rev_diff_weight,sizeof(unsigned int)*E) );


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
        }
   
    }

    getCx<<<1,1>>>(endNode,D_dest_cost);
    gpuErrchk( cudaMemcpy(H_dest_cost,D_dest_cost, sizeof(int),cudaMemcpyDeviceToHost) );
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
        Path.push_back(p);
        printf("%d\n",p);
    }
    else{
        printf("not found\n");
    }

    //reverse the path to get from source to end
    reverse(Path.begin(),Path.end());

    
    ///////////////////////////////////////////////
    // A star complete //

    FILE* fdiff = fopen("Updates.txt","r");
    int line;
    while(fscanf(fdiff,"%d\n",&line)!=EOF){
        unordered_map<unsigned int,Node*> Graph;
        unordered_map<unsigned int,Node*> rev_Graph;

        bool flag_do_a_star = false;
        
        int insertEdge=0, delEdge=0;
        int delEdgesV_size = 0;                    //v whose cost can change due to deletion
        for(int i=0;i<line;i++){
            int flag;
            int u,v;
            unsigned int w;
            fscanf(fdiff,"%d %d %d %u\n",&flag,&u,&v,&w);
            if(flag==1){
                insertDiff(Graph,u,v,w);
                insertDiff(rev_Graph,v,u,w);
                insertEdge++;
            }
            else if(flag==0){
                //check id del edges in optimal path.
                check_del_path(u,v,Path,flag_do_a_star);
                removeDelEdges(u,v,H_offset,H_edges,N,E,H_rev_offset,H_rev_edges);
                //add to list only if its cost changes due to this deletion
                if(H_parent[v]==u){
                    H_delEdgesV[delEdgesV_size]=v;
                    delEdgesV_size++;
                }
                
                delEdge++;
            }
            
        }

        if(1)
            printf("[INFO] insertion:%d, deletion:%d\n",insertEdge,delEdge);

        //reset offset to 0 ..ie no nodes
        memset(H_diff_offset,0,sizeof(int)*N);
        memset(H_rev_diff_offset,0,sizeof(int)*N);

        createDiffGraph(N,Graph,H_diff_offset,H_diff_edges,H_diff_weight);
        createDiffGraph(N,rev_Graph,H_rev_diff_offset,H_rev_diff_edges,H_rev_diff_weight);
        
        //TODO free the graphs

        
        //deleted edges
        gpuErrchk ( cudaMemcpy(D_edges,H_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_edges,H_rev_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_delEdgesV,H_delEdgesV,sizeof(int)*E,cudaMemcpyHostToDevice) );

        //diff graph
        gpuErrchk ( cudaMemcpy(D_diff_edges,H_diff_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_diff_offset,H_diff_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_diff_weight,H_diff_weight,sizeof(unsigned int)*E,cudaMemcpyHostToDevice) );

        //rev diff graph
        gpuErrchk ( cudaMemcpy(D_rev_diff_edges,H_rev_diff_edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_diff_offset,H_rev_diff_offset,sizeof(int)*N,cudaMemcpyHostToDevice) );
        gpuErrchk ( cudaMemcpy(D_rev_diff_weight,H_rev_diff_weight,sizeof(unsigned int)*E,cudaMemcpyHostToDevice) );

        //reset D_nV flag
        gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

        //add del
        if(delEdgesV_size>0){
            if(DEBUG)
                printf("[INFO] Starting computing cost for deletions\n");

            int numBlocks_del = ( delEdgesV_size + numThreads -1)/numThreads;
            propogateDel<<<numBlocks_del,numThreads>>>(D_delEdgesV,delEdgesV_size,D_rev_offset,D_rev_edges,D_rev_weight,N,E,
                D_hx,D_parent,D_lock,D_nVFlag);
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();
        
        }

        if(DEBUG)
            printf("[INFO] starting computing cost for inserions\n");
        
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

        //copy back
        gpuErrchk( cudaMemcpy(H_nV_size,D_nV_size, sizeof(int),cudaMemcpyDeviceToHost) );

        //reset nV flags
        gpuErrchk( cudaMemcpy(D_nVFlag,H_nVFlag,sizeof(int)*N,cudaMemcpyHostToDevice) );

        if(DEBUG)
            printf("[INFO] startig propogation\n");

        while(*H_nV_size > 0){

            numBlocks = (*H_nV_size+numThreads-1)/numThreads;

            propogate<<<numBlocks,numThreads>>>(D_nV,D_nV_size,D_offset,D_edges,D_weight,D_hx,
                N,E,D_lock,D_parent,D_nVFlag,
                D_diff_offset,D_diff_edges,D_diff_weight,insertEdge,
                D_rev_offset,D_rev_edges,D_rev_weight,
                D_rev_diff_offset,D_rev_diff_edges,D_rev_diff_weight);
            
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


        if(DEBUG)
            printf("[INFO] updating priority queue\n");

        //propogate complete do normal A*
        numBlocks = (K+numThreads-1)/numThreads;

        //update PQ after propogate
        keepHeapPQ<<<numBlocks,numThreads>>>(D_PQ_size,N,K);
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();

        //check if there is node cost in PQ less than dest
        *H_flagEnd = 1;
        gpuErrchk( cudaMemcpy(D_flagEnd,H_flagEnd,sizeof(int),cudaMemcpyHostToDevice) );
    
        checkMIN<<< numBlocks,numThreads >>>(D_PQ_size,D_flagEnd,endNode,N,K);
        
        gpuErrchk( cudaPeekAtLastError() );

        gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );

        //here flag end represents from above that there is a node with cost lesser     
        if(*H_flagEnd==0 && flag_do_a_star){

            printf("[INFO] doing a* after propogation\n");
            
            insertDest<<<1,1>>>(D_PQ_size,endNode,D_openList);
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();

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
                extractMin<<<numBlocks,numThreads>>>(D_PQ_size, D_expandNodes,D_expandNodes_size,D_openList,N,K);
                
                gpuErrchk(cudaPeekAtLastError() );
    
                cudaDeviceSynchronize();
    
                
                A_star_expand<<<numBlocks,numThreads>>>(D_offset,D_edges,D_weight,D_hx,D_parent,
                    D_expandNodes,D_expandNodes_size, D_lock ,D_flagfound,D_openList,
                    N,E,K,endNode,D_nVFlag,D_PQ_size,
                    true,D_diff_offset,D_diff_edges,D_diff_offset,insertEdge);
                
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
                    
                    gpuErrchk( cudaPeekAtLastError() );
                    cudaDeviceSynchronize();
                    gpuErrchk( cudaMemcpy(H_flagEnd,D_flagEnd, sizeof(int),cudaMemcpyDeviceToHost) );
                // printf("\ninside MIN\n");
                }
        
            }
        }
    
        getCx<<<1,1>>>(endNode,D_dest_cost);
        gpuErrchk( cudaMemcpy(H_parent,D_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );


        // found or not found based on Cx
        gpuErrchk( cudaMemcpy(H_dest_cost,D_dest_cost, sizeof(int),cudaMemcpyDeviceToHost) );

        printf("[OUT] Updated Cost: %d\n",*H_dest_cost);
        printf("[OUT] PATH(in reverse): ");
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

void removeDelEdges(int u,int v,int* offset,int* edges,int N,int E,int* rev_offset,int* rev_edges){
    int start = offset[u];
    int end =  E;
    if(u!=N-1)
        end = offset[u+1];
    while(start<end){
        if( v == edges[start]){
            edges[start]=-1;
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
            break;
        }
        start++;
    }

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