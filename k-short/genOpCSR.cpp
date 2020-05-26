#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unordered_map>

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

void insertG(unordered_map< unsigned int, Node*>& Graph,int a,int b,int c);
void printGraphCSR(int N,int E,unordered_map<unsigned int,Node*>& Graph, string filename);

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

    unordered_map<unsigned int,Node*> Graph;

    for(int i=0;i<N;i++){
        int start = H_offset[i];
        int end = E;
        if(i != N-1){
            end = H_offset[i+1];
        }

        while (start < end)
        {
            int child = H_edges[start];
            int weight = H_weight[start];

            insertG(Graph,child,i,weight);
            start++;
        }
        

    }

    printGraphCSR(N,E,Graph,"graph_op.txt");

    return 0;
}

void insertG(unordered_map< unsigned int, Node*>& Graph,int a,int b,int c){
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


void printGraphCSR(int N,int E,unordered_map<unsigned int,Node*>& Graph, string filename){
    FILE* fptr = fopen(filename.c_str(),"w");

    fprintf(fptr,"%d %d\n",N,E);
    
    vector<unsigned int> off;
    off.push_back(0);
    int k =0;
    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->Edges.size();j++){
                fprintf(fptr,"%d ",n->Edges[j]->val);
                k++;
            }
            off.push_back(k);
        }
        else{
            off.push_back(k);
        }
    }

    fprintf(fptr,"\n");

    for(int j=0;j<off.size()-1;j++){
        fprintf(fptr,"%d ",off[j]);
    }

    fprintf(fptr,"\n");
    

    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->weights.size();j++){
                fprintf(fptr,"%d ",n->weights[j]);
                k++;
            }
        }
    }
    
    fprintf(fptr,"\n");
    fclose(fptr);

}
