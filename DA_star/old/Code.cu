#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(){
    int N,E;
    scanf("%d %d\n",&N,&E);
    int* H_Edges = (int*)malloc(sizeof(int)*E);
    int* H_Offset = (int*)malloc(sizeof(int)*N);
    unsigned int* H_Weight = (unsigned int*)malloc(sizeof(unsigned int)*E);

    for(int i=0;i<E;i++){
        scanf("%d",&H_Edges[i]);
    }

    for(int i=0;i<N;i++){
        scanf("%d",&H_Offset[i]);
    }

    for(int i=0;i<E;i++){
        scanf("%u",&H_Weight[i]);
    }

}
