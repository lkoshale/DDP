#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <cstdlib>
#include <unordered_map>
#include <vector>

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

    void remove_edge(Node* v){
        vector<Node*>::iterator it = Edges.begin();
        while(it!= Edges.end()){
            if((*it)->val == v->val){
                Edges.erase(it);
                break;
            }

            it++;
        }
            
    }

};


void printGraphCSR(int N,int E, vector<Node*>& Graph, string filename);

void print_maze(vector<bool>& obstacle_flag,int N,int M, FILE* fptr);

int main()
{
    int N;
    int M;
    cin>>N>>M;

    float prob;
    cin>>prob;

    int num_nodes = N*M;

    vector<Node*> Graph;

    vector< int > obstacle;

    vector<bool> obstacle_flag;

    for(int i=0;i<num_nodes;i++){
        Node* temp = new Node(i);
        Graph.push_back(temp);
    }


    for(int i=0;i<num_nodes;i++){
        int r = rand()%100;
        if(prob*100 > r ){
            obstacle.push_back(i);
            obstacle_flag.push_back(true);
        }
        else{
            obstacle_flag.push_back(false);
        }
        
    }


    int E = 0;

    for(int i=0;i<N;i++){
        
        for(int j=0;j<M;j++){

            Node* node = Graph[i*M + j];

            if(obstacle_flag[i*M+j])
                continue;

            if(i > 0  && (!obstacle_flag[(i-1)*M+j]) ){
                node->addEdge(Graph[(i-1)*M+j],1);
                E++;
            }

            if(i < N-1 && (!obstacle_flag[(i+1)*M+j]) ){
                node->addEdge(Graph[(i+1)*M+j],1);
                E++;
            }

            if(j < M - 1 && (!obstacle_flag[i*M+(j+1)] ) ){
                node->addEdge(Graph[i*M+(j+1)],1);
                E++;
            }
            
            if(j > 0 && (!obstacle_flag[i*M+(j-1)] ) ){
                node->addEdge(Graph[i*M+(j-1)],1);
                E++;
            }
        }
    }

    printGraphCSR(num_nodes,E,Graph,"graph.txt");

    FILE* fviz = fopen("viz_graph.txt","w");
    print_maze(obstacle_flag,N,M,fviz);

    

    FILE* fhx = fopen("Hx.txt","w");
    int e_i = N-1;
    int e_j = M-1;
    for(int i=0;i<num_nodes;i++){
        int t_i = i / M;
        int t_j = i % M;
        int dist = (e_i - t_i)*(e_i - t_i) + (e_j-t_j)*(e_j-t_j);
        fprintf(fhx,"%d ",dist);
    }
    fclose(fhx);
    

    int counter = 10;
    cin>>counter;

    FILE* fupd = fopen("Updates.txt","w");

    for(int cn=0;cn<counter;cn++){

        //or move obstacle
        int c_size = obstacle.size()/counter;

        int n_E=0;
        string str = "";
        for(int i=0;i<c_size;i++){
            int n = obstacle[obstacle.size()-1];
            obstacle.pop_back();
            obstacle_flag[n] = false;
            
            Node* node = Graph[n];
            int c_i = n / M;
            int c_j = n % M;


            if(c_i > 0  ){
                node->addEdge(Graph[(c_i-1)*M+c_j],1);
                Graph[(c_i-1)*M+c_j]->addEdge(node,1);
                
                str+="1 "+to_string(n)+" "+to_string((c_i-1)*M+c_j)+" 1\n";
                str+="1 "+to_string((c_i-1)*M+c_j)+" "+to_string(n)+" 1\n";
            
                n_E+=2;
            }

            if(c_i < N-1  ){
                node->addEdge(Graph[(c_i+1)*M+c_j],1);
                Graph[(c_i+1)*M+c_j]->addEdge(node,1);

                str+="1 "+to_string(n)+" "+to_string((c_i+1)*M+c_j)+" 1\n";
                str+="1 "+to_string((c_i+1)*M+c_j)+" "+to_string(n)+" 1\n";
                n_E+=2;
            }

            if(c_j < M - 1  ){
                node->addEdge(Graph[c_i*M+(c_j+1)],1);
                Graph[c_i*M+(c_j+1)]->addEdge(node,1);

                str+="1 "+to_string(n)+" "+to_string(c_i*M+(c_j+1))+" 1\n";
                str+="1 "+to_string(c_i*M+(c_j+1))+" "+to_string(n)+" 1\n";
                n_E+=2;
            }
            
            if(c_j > 0 ){
                node->addEdge(Graph[c_i*M+(c_j-1)],1);
                Graph[c_i*M+(c_j-1)]->addEdge(node,1);
                str+="1 "+to_string(n)+" "+to_string(c_i*M+(c_j-1))+" 1\n";
                str+="1 "+to_string(c_i*M+(c_j-1))+" "+to_string(n)+" 1\n";

                n_E+=2;
            }
        }


        fprintf(fupd,"%d\n",n_E);
        fprintf(fupd,"%s",str.c_str());



    }



}



void printGraphCSR(int N,int E,vector<Node*>& Graph, string filename){
    FILE* fptr = fopen(filename.c_str(),"w");

    fprintf(fptr,"%d %d\n",N,E);
    
    vector<unsigned int> off;
    off.push_back(0);
    int k =0;
    for(int i=0;i<N;i++){

        Node* n = Graph[i];
        for(int j=0;j<n->Edges.size();j++){
            if(n->Edges[j]!=NULL){
                fprintf(fptr,"%d ",n->Edges[j]->val);
                k++;
            }
        }
        off.push_back(k);
       
    }

    fprintf(fptr,"\n");

    for(int j=0;j<off.size()-1;j++){
        fprintf(fptr,"%d ",off[j]);
    }

    fprintf(fptr,"\n");
    

    for(int i=0;i<N;i++){
        
        Node* n = Graph[i];
        for(int j=0;j<n->weights.size();j++){
            if(n->Edges[j]!=NULL){
                fprintf(fptr,"%d ",n->weights[j]);
            }
        }
    
    }
    
    fprintf(fptr,"\n");
    fclose(fptr);

}


void print_maze(vector<bool>& obstacle_flag,int N,int M, FILE* fptr)
{
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            if(obstacle_flag[i*M+j])
                fprintf(fptr,".");
            else
                fprintf(fptr,"*");
        }

        fprintf(fptr,"\n");
    }

}