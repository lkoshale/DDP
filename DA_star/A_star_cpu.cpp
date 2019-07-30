#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <climits>
#include <cstdio>

using namespace std;

void insertPQ(int* PQ,int* size,int node,int* Cx){
    if(*size==0){
        PQ[*size]=node;
        (*size)++;
    }
    else{
        PQ[*size]=node;
        (*size)++;

        int i = *size - 1;
        while(i>0){
            if(Cx[PQ[(i-1)/2]] > Cx[PQ[i]] ){
                int swap = PQ[i];
                PQ[i] = PQ[(i-1)/2];
                PQ[(i-1)/2] = swap;
                i = (i-1)/2;
            }
            else
                break;
        }
    }
}

int getMin(int* PQ,int* size,int* Cx){
    if(*size==0)
        return -1;
    
    int min = PQ[0];

    //swap with last
    PQ[0]= PQ[*size -1];
    (*size)--;
    int i=0;
    while(2*i+1< *size){
        if(2*i+2 >= *size){
            if(Cx[PQ[2*i+1]] < Cx[PQ[i]]){
                int swap = PQ[2*i+1];
                PQ[2*i+1]=PQ[i];
                PQ[i] = swap;
                i = 2*i +1;
            }
            else 
                break;
        }
        else{
            if( Cx[PQ[i]] > Cx[PQ[2*i+1]] && Cx[PQ[2*i+1]]<= Cx[PQ[2*i+2]]){
                int swap = PQ[2*i+1];
                PQ[2*i+1]=PQ[i];
                PQ[i] = swap;
                i = 2*i +1;
            }
            else if(Cx[PQ[i]] > Cx[PQ[2*i+2]] && Cx[PQ[2*i+2]]<= Cx[PQ[2*i+1]]){
                int swap = PQ[2*i+2];
                PQ[2*i+2]=PQ[i];
                PQ[i] = swap;
                i = 2*i +2;
            }
            else
                break;
        }

    }

    return min;
}

//update happens when Cx is decreased
void updatePQ(int* PQ,int* size,int node,int* Cx){
    int index=-1;
    for(int i=0;i<*size;i++){
        if(PQ[i]==node)
            index = i;
    }

    if(index>0){
        int i = index;
        while(i>0){
            if(Cx[PQ[(i-1)/2]] > Cx[PQ[i]] ){
                int swap = PQ[i];
                PQ[i] = PQ[(i-1)/2];
                PQ[(i-1)/2] = swap;
                i = (i-1)/2;
            }
            else
                break;
        }
    }
}

//void update(int*)
int main(){
    // set start node same as destination
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
    int* H_openList  = (int*)malloc(sizeof(int)*N);


    int* PQ = (int*)malloc(sizeof(int)*N);
    int* PQ_size = (int*)malloc(sizeof(int));
    *PQ_size = 0;


    memset(H_parent,-1,sizeof(int)*N);
    memset(H_openList,-1,sizeof(int)*N);

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

    //A*
    H_cx[startNode]=H_hx[startNode];
    insertPQ(PQ,PQ_size,startNode,H_cx);
    H_openList[startNode]=1;
    bool flag_found = false;
    int curr_node ;
    while( *PQ_size >0 && !flag_found){
        curr_node = getMin(PQ,PQ_size,H_cx);
        //remove from openList
        H_openList[curr_node]=-1;

        //got dest
        if(curr_node==endNode){
            flag_found = true;
            break;
        }

        //iterate 
        int begin = H_offset[curr_node];
        int end = E-1;
        if(curr_node!=N-1)
            end = H_offset[curr_node+1];
        
        for(int j=begin;j<end;j++){
            int child = H_edges[j];
            if(H_cx[child] > H_cx[curr_node]-H_hx[curr_node]+H_weight[j]+H_hx[child]){
                H_cx[child] = H_cx[curr_node] - H_hx[curr_node] + H_weight[j]+H_hx[child];
                if(H_openList[child]==1){
                    updatePQ(PQ,PQ_size,child,H_cx);
                }
                else{
                    insertPQ(PQ,PQ_size,child,H_cx);
                    H_openList[child]=1;
                }
                H_parent[child]=curr_node;
            }
        }
    }

    if(flag_found){
        printf("cost %d\n",H_cx[endNode]);

        int c_node = endNode;
        printf("%d ",c_node);
        while(H_parent[c_node]!=-1){
            printf("%d ",H_parent[c_node]);
            c_node = H_parent[c_node];
        }
        printf("\n");
    }
    else{
        printf("unreachable !!\n");
    }    
}