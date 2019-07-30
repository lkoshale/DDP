#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
 
#define  MAX_EDGE  100000000
 

__global__ void BFS(int* off,int* edge,int* current,int* size,int N,int E,int* c_arr,int* c_size,int* dist){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < *size){
        // printf("inside kernel %d %d\n",id,*size);
        int node = current[id];
        //extend this node
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start<end){
            // add in list atomically in c_arr
            int child = edge[start];
            // printf("c %d \n",child);
            if ( dist[child] < 0){
                dist[child] = dist[node] + 1;
                int index = atomicAdd(c_size,1);
                c_arr[index]= child;
            }
            start++;  
        }

        // printf("s %d\n",*c_size);
    }
}


int main(){

    // set start node same as destination
    int startNode;
    scanf("%d",&startNode);

    FILE* fgraph = fopen("graph.txt","r");
    int N,edge_size;
    fscanf(fgraph,"%d %d\n",&N,&edge_size);


    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*edge_size);


    for(int i=0;i<edge_size;i++){
        fscanf(fgraph,"%d",&H_edges[i]);
    }

    for(int i=0;i<N;i++){
        fscanf(fgraph,"%d",&H_offset[i]);
    }

    
    
    fclose(fgraph);
    printf("completed input\n");

    int* H_current_node = (int*)malloc(sizeof(int)*edge_size);
    H_current_node[0]=startNode;
    int* H_c_size = (int*)malloc(sizeof(int));
    *H_c_size = 1;
    int* H_visited = (int*)malloc(sizeof(int)*N);
    memset(H_visited,-1,sizeof(int)*N);
    H_visited[startNode]=0;

    int* a0 = (int*)malloc(sizeof(int));
    *a0=0;

    int* a1 = (int*)malloc(sizeof(int));
    *a1=1;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    int* D_offset;
    int* D_edges;
    int* D_visited;

    int* D_current_node1;
    int* D_c_size1;
    int* D_current_node2;
    int* D_c_size2;

    cudaMalloc(&D_offset,sizeof(int)*N);
    cudaMalloc(&D_visited,sizeof(int)*N);
    cudaMalloc(&D_edges,sizeof(int)*edge_size);
    cudaMalloc(&D_current_node1,sizeof(int)*edge_size);
    cudaMalloc(&D_c_size1,sizeof(int));
    cudaMalloc(&D_current_node2,sizeof(int)*edge_size);
    cudaMalloc(&D_c_size2,sizeof(int));


    cudaMemcpy(D_offset,H_offset,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_edges,H_edges,sizeof(int)*edge_size,cudaMemcpyHostToDevice);
    cudaMemcpy(D_current_node1,H_current_node,sizeof(int)*edge_size,cudaMemcpyHostToDevice);
    cudaMemcpy(D_visited,H_visited,sizeof(int)*N,cudaMemcpyHostToDevice);
    cudaMemcpy(D_c_size1,a1,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(D_c_size2,a0,sizeof(int),cudaMemcpyHostToDevice);



    int i=1;
    cudaEventRecord(start);
    while(*H_c_size>0){
        int numThreads = 512;
        int numBlocks = (*H_c_size+numThreads-1)/numThreads;


        if(i%2==1){
            //use array 1
            BFS<<<numBlocks,numThreads>>>(D_offset,D_edges,D_current_node1,D_c_size1,N,edge_size,D_current_node2,D_c_size2,D_visited);

            cudaMemcpy(H_c_size,D_c_size2, sizeof(int),cudaMemcpyDeviceToHost);
            // reset the index
            cudaMemcpy(D_c_size1,a0,sizeof(int),cudaMemcpyHostToDevice);
          
        }
        else{
            //use array 2
            BFS<<<numBlocks,numThreads>>>(D_offset,D_edges,D_current_node2,D_c_size2,N,edge_size,D_current_node1,D_c_size1,D_visited);
            
            cudaMemcpy(H_c_size,D_c_size1, sizeof(int),cudaMemcpyDeviceToHost);
            //reset index
            cudaMemcpy(D_c_size2,a0,sizeof(int),cudaMemcpyHostToDevice);

        }
        i++;
        
    }
    cudaEventRecord(stop);
    
    cudaMemcpy(H_visited,D_visited, sizeof(int)*N,cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    int max_level = 0;
    for(int i=0;i<N;i++){
        if(H_visited[i]>max_level){
            max_level = H_visited[i];
        }
        // printf("%d, %d\n",i,H_visited[i]);
    }
    
    printf("max-level: %d\n",max_level);
    printf("time: %f\n",milliseconds);

    FILE* fptr = fopen("Hx.txt","w");
    int count_unreachable = 0;
    int farthest_val = 0;
    int farthest_node = 0;
    for(int i=0;i<N;i++){
        if(H_visited[i]==-1)
            count_unreachable++;
        if(H_visited[i]>=farthest_val){
            farthest_val = H_visited[i];
            farthest_node = i;
        }

        fprintf(fptr,"%d ",H_visited[i]);
    }
    fclose(fptr);

    printf("unreachable %d , %f\n",count_unreachable,count_unreachable/(N*1.0));
    printf("farthest node: %d  val: %d\n",farthest_node,farthest_val);
    return 0;
}
