#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stack>


using namespace std;

class Node{
public:
    unsigned int val;
    vector<unsigned int> weights;
    vector<Node*> Edges;
    vector<Node*> backEdge;
    vector<unsigned int> backWeight;

    Node(int val){
        this->val = val;
    }

    void addEdge(Node* v,unsigned int w){
        this->Edges.push_back(v);
        this->weights.push_back(w);
    }

    void addBackEdge(Node* v, unsigned int w){
        this->backEdge.push_back(v);
        this->backWeight.push_back(w);
    }

};

void genUpdates(FILE* fptr,int N, int& E, unordered_map<unsigned int, Node*>& Graph,int K);

void printGraphCSR(int N,int E,unordered_map<unsigned int,Node*>& Graph, string filename);

void printBackCSR(int N, int E, unordered_map<unsigned int,Node*>& Graph,string filename);




int main(){

    srand(42);

    int K = 4;

    unordered_map<unsigned int,Node*> Graph;

    unsigned int N = 0;
    unsigned int a,b,c;

    int E = 0;

    while( scanf("%u %u %u\n",&a,&b,&c)!=EOF){
        // c = rand()%100;
        
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(a);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            unordered_map<unsigned int,Node*>:: iterator it;
            it = Graph.find(b);
            if(it!=Graph.end()){
                Node* v = it->second;
                n->addEdge(v,c);
                v->addBackEdge(n,c);
            }
            else{
                Node* v = new Node(b);
                n->addEdge(v,c);
                v->addBackEdge(n,c);
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
                v->addBackEdge(n,c);
            }
            else{
                Node* v = new Node(b);
                n->addEdge(v,c);
                v->addBackEdge(n,c);
                Graph.insert(pair<unsigned int,Node*>(b,v));
            }

        }

        if( a > N ) N=a;
        if( b > N ) N=b;

        E++;
    }
    N++;
    FILE* fptr = fopen("Updates.txt","w");
    for(int i=0;i<K;i++)
        genUpdates(fptr,N,E,Graph,K);
    
    fclose(fptr);

    printGraphCSR(N,E,Graph,"graph.txt");

   // printBackCSR(N,E,Graph,"graph_op.txt");



    return 0;
}

void genUpdates(FILE* fptr,int N, int& E, unordered_map<unsigned int, Node*>& Graph,int K){
    int remove = E/(2*K);
    int count = 0;
    string strEdges = "";
    for(int i=0;i<remove/2;i++){
        int a = rand()%N;
        unordered_map<unsigned int, Node*>::iterator itr;
        itr = Graph.find(a);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            if(n->Edges.size() > 0){
                int b = rand()%(n->Edges.size());
                Node* m = n->Edges[b];
                if(m!=NULL){
                    strEdges += "1 "+to_string(n->val)+" "+to_string(m->val)+" "+to_string(n->weights[b])+"\n";
                    
                    n->Edges[b] =  NULL;
                    count++;
                }
                
            }
        }
    }


    int delcount = 0;
    /*
    for(int i=0;i<remove/2;i++){
        int a = rand()%N;
        unordered_map<unsigned int, Node*>::iterator itr;
        itr = Graph.find(a);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            if(n->Edges.size() > 0){
                int b = rand()%(n->Edges.size());
                Node* m = n->Edges[b];
                if(m!=NULL){
                    strEdges += "0 "+to_string(n->val)+" "+to_string(m->val)+" "+to_string(n->weights[b])+"\n";
                    delcount++;
                }
                
            }
        }
    }*/



    E = E-count;
    fprintf(fptr,"%d\n",count+delcount);
    fprintf(fptr,"%s",strEdges.c_str());

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
    

    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->weights.size();j++){
                if(n->Edges[j]!=NULL){
                    fprintf(fptr,"%d ",n->weights[j]);
                }
            }
        }
    }
    
    fprintf(fptr,"\n");
    fclose(fptr);

}

void printBackCSR(int N, int E, unordered_map<unsigned int,Node*>& Graph,string filename){
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
            for(int j=0;j<n->backEdge.size();j++){
                fprintf(fptr,"%d ",n->backEdge[j]->val);
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
            for(int j=0;j<n->backWeight.size();j++){
                fprintf(fptr,"%d ",n->backWeight[j]);
                k++;
            }
        }
    }
    
    fprintf(fptr,"\n");

    fclose(fptr);
}