#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <queue>

#define  MAX_EDGE  100000000

using namespace std;

// class Compare{
//     int* 
// public:

//     bool operator()(int a,int b){
//         return true;
//     }
// };


int main(){
    int N;
    scanf("%d\n",&N);
    int startNode;
    scanf("%d\n",&startNode);
    int endNode;
    scanf("%d\n",&endNode);

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
    memset(H_cx,-1,sizeof(int)*N);
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
    
    for(int i=0;i<N;i++){
        printf("%d ",H_offset[i]);
    }

    printf("\n");
    
    for(int i=0;i<edge_size;i++){
        printf("%d ",H_edges[i]);
    }

    printf("\n");
    
    for(int i=0;i<edge_size;i++){
        printf("%d ",H_weights[i]);
    }
    


    return 0;

}