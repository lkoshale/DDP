#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>

#define  MAX_EDGE  10000000

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


void insertDiff(unordered_map< unsigned int, Node*>& Graph,int a,int b,unsigned int c);

void createDiffGraph(int N,unordered_map<unsigned int,Node*>& Graph,
    int* diffOff,int* diffEdges,unsigned int* diffWeight );

void removeDelEdges(int u,int v,int* offset,int* edges,int N,int E);

void mergeDiff(int* offset,int* edges,unsigned int* weight,int N,int& E,
            int* diff_offset, int* diff_edges,unsigned int* diff_weight,int insert_size,int del_size,
            int* mOffset,int* mEdges,unsigned int* mWeight);

int main(){
    FILE* fgraph = fopen("graph.txt","r");
    int N,E;
    fscanf(fgraph,"%d %d\n",&N,&E);

    int* H_offset = (int*)malloc(sizeof(int)*N);
    int* H_edges  = (int*)malloc(sizeof(int)*E);
    unsigned int* H_weight = (unsigned int*)malloc(sizeof(unsigned int)*E);

    
    for(int i=0;i<E;i++){
        fscanf(fgraph,"%d",&H_edges[i]);
    }

    for(int i=0;i<N;i++){
        fscanf(fgraph,"%d",&H_offset[i]);
    }

    for(int i=0;i<E;i++){
        fscanf(fgraph,"%u",&H_weight[i]);
    }

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

                removeDelEdges(u,v,H_offset,H_edges,N,E);
                
                delEdge++;
            }
            
        }

        int* H_diff_edges = (int*)malloc(sizeof(int)*insertEdge);
        int* H_diff_offset = (int*)malloc(sizeof(int)*N);
        unsigned int* H_diff_weight = (unsigned int*)malloc(sizeof(unsigned int)*insertEdge);
         //reset offset to 0 ..ie no nodes
        memset(H_diff_offset,0,sizeof(int)*N);

        if(1)
            printf("[INFO] insertion:%d, deletion:%d, delaff:%d\n",insertEdge,delEdge,delEdgesV_size);

       

        createDiffGraph(N,Graph,H_diff_offset,H_diff_edges,H_diff_weight);
        for(int i=0;i<insertEdge;i++){
            printf("%d ",H_diff_edges[i]);
        }
        printf("\n");
        for(int i=0;i<N;i++){
            printf("%d ",H_diff_offset[i]);
        }
        printf("\n");
        for(int i=0;i<insertEdge;i++){
            printf("%d ",H_diff_weight[i]);
        }
        printf("\n\n");

        int* mOffset,*mEdges;
        unsigned int* mWeight;
        
        int mE = E + insertEdge - delEdge;

        mOffset = (int*)malloc(sizeof(int)*N);
        mEdges =  (int*)malloc(sizeof(int)*mE);
        mWeight = (unsigned int*)malloc(sizeof(unsigned int)*mE);

        mergeDiff(H_offset,H_edges,H_weight,N,E,
            H_diff_offset,H_diff_edges,H_diff_weight,insertEdge,delEdge,
            mOffset,mEdges,mWeight);

        
        free(H_offset);
        free(H_edges);
        free(H_weight);
        free(H_diff_offset);
        free(H_diff_edges);
        free(H_diff_weight);

        H_offset = mOffset;
        H_edges = mEdges;
        H_weight = mWeight;
        E = mE;

        

    }

        for(int i=0;i<E;i++){
            printf("%d ",H_edges[i]);
        }
        printf("\n");
        for(int i=0;i<N;i++){
            printf("%d ",H_offset[i]);
        }
        printf("\n");
        for(int i=0;i<E;i++){
            printf("%d ",H_weight[i]);
        }
        printf("\n");

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

void removeDelEdges(int u,int v,int* offset,int* edges,int N,int E){
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


}


void mergeDiff(int* offset,int* edges,unsigned int* weight,int N,int& E,
            int* diff_offset, int* diff_edges,unsigned int* diff_weight,int insert_size,int del_size,
            int* mOffset,int* mEdges,unsigned int* mWeight){

  
    int mE = E + insert_size - del_size;
    mOffset[0] = 0;
    int edegOffset= 0;
    for(int i=0;i<N;i++){

        int start = offset[i];
        int end = E;
        if(i!=N-1)
            end = offset[i+1];
        int count = 0;
        while(start<end){
            int child = edges[start];
            if(child!=-1){
                mEdges[edegOffset+count] = child;
                mWeight[edegOffset+count] = weight[start];
                count++;
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
                mEdges[edegOffset+count] = child;
                mWeight[edegOffset+count]= diff_weight[start];
                count++;
            }
            start++;
        }

        edegOffset+=count;
        if(i!=N-1)
            mOffset[i+1]=edegOffset;
    }

    //change E


    for(int i=0;i<mE;i++){
        printf("%d ",mEdges[i]);
    }
    printf("\n");
    for(int i=0;i<N;i++){
        printf("%d ",mOffset[i]);
    }
    printf("\n");
    for(int i=0;i<mE;i++){
        printf("%d ",mWeight[i]);
    }
    printf("\n");

}