#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>

using namespace std;

class Node{
public:
    unsigned int val;
    vector<Node*> Edges;

    Node(int val){
        this->val = val;
    }

    void addEdge(Node* v){
        this->Edges.push_back(v);
    }

};


void printGraphCSR(int N,int E,unordered_map<unsigned int,Node*>& Graph, string filename);

int main(){

    vector<float> x;
    vector<float> y;

    FILE* fcord = fopen("Cord.txt","r");

    float cx,cy;
    while(fscanf(fcord,"%f %f \n",&cx,&cy)!=EOF){
        x.push_back(cx);
        y.push_back(cy);
    }


    //create graph
    int range = 15;
    int N = x.size();
    int E = 0;

    unordered_map<unsigned int,Node*> Graph;
    for(unsigned int i=0;i<N;i++){
        Node* t = new Node(i);
        Graph[i] = t;
    }

    for(unsigned int i=0;i<N-1;i++){
        float cx,cy;
        cx= x[i];
        cy= y[i];
        Node* u = Graph[i];
        for(unsigned int j=0;j<N;j++){
            if(j==i)
                continue;
            float dist_sq = (cx-x[j])*(cx-x[j]) + (cy-y[j])*(cy-y[j]);
            if(dist_sq < range*range){
                Node* v = Graph[j];
                u->addEdge(v);
                E++;
            }
        }

    }


    printGraphCSR(N,E,Graph,"graph.txt");

    return 0;
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
                if(n->Edges[j]!=NULL){
                    fprintf(fptr,"%d ",n->Edges[j]->val);
                    k++;
                }
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
    fclose(fptr);

}