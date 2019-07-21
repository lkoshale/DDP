#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <unordered_map>
#include<time.h>

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

int main(){

    srand(42);

    unordered_map<unsigned int,Node*> Graph;

    unsigned int N = 0;
    unsigned int a,b,c;

    int E = 0;

    while( scanf("%u %u\n",&a,&b,&c)!=EOF){
        c = rand()%100;
        
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

        if( a > N ) N=a;
        if( b > N ) N=b;

        E++;
    }
    N++;

    cout<<N << " " << E <<"\n";

    vector<unsigned int> off;

    off.push_back(0);
    int k =0;
    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->Edges.size();j++){
                cout<<n->Edges[j]->val<<" ";
                k++;
            }
            off.push_back(k);
        }
        else{
            off.push_back(k);
        }
    }

    cout<<"\n";
    for(int j=0;j<off.size()-1;j++){
        cout<<off[j]<<" ";
    }

    cout<<"\n";

    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->weights.size();j++){
                cout<<n->weights[j]<<" ";
                k++;
            }
        }
    }
    cout<<"\n";


    return 0;
}