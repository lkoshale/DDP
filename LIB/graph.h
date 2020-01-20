
#ifndef GRAPH_H
#define GRAPH_H


#include <unordered_map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "graph_builder.h"

template <class DataType , class FuncType = Builder<DataType> > 
class Graph 
{
    private:
        unsigned int N;         //nodes from 0 to N-1
        unsigned int E;
        int* offsets;
        int* edges;
        DataType* weights;

        std::unordered_map<long,int> vertex_map;     //if vertex doesn't follow 0-1 naming

    public: 
        Graph();

        void __alloc__(int n, int e);

        void build_graph( FuncType function );

        int get_num_nodes();
        int get_num_edges();

        int* get_offsets();
        int* get_edges();
        DataType* get_weight();

        void read_graph(std::string filename);

        void free();

        
};



/*
* Implementation Here as 
* templates cant have diff files for implementation
* (or move to another file .tpp nand include it here)
*/

template <class DataType, class FuncType >
Graph<DataType, FuncType>::Graph()
{
    N = 0;
    E = 0;
}

template <class DataType, class FuncType >
void Graph<DataType, FuncType> :: __alloc__(int n, int e){
    N = n;
    E = e;
    this->offsets = (int*)malloc(sizeof(int)*N);
    this->edges = (int*)malloc(sizeof(int)*E);
    this->weights = (DataType*)malloc(sizeof(DataType)*E);
}


template <class DataType, class FuncType >
void Graph<DataType, FuncType> :: build_graph(FuncType builder )
{
    builder(this->offsets,this->edges,this->weights);
}

template <class DataType, class FuncType >
int Graph<DataType, FuncType> :: get_num_nodes(){ return N; }

template <class DataType, class FuncType >
int Graph<DataType, FuncType> :: get_num_edges(){ return E; }

template <class DataType, class FuncType >
int* Graph<DataType, FuncType> :: get_offsets(){ return this->offsets; }

template <class DataType, class FuncType >
int*  Graph<DataType, FuncType> :: get_edges(){ return this->edges; }

template <class DataType, class FuncType >
DataType*  Graph<DataType, FuncType> :: get_weight(){ return this->weights; }


template <class DataType, class FuncType >
void Graph<DataType, FuncType> :: read_graph(std::string filename){
    std:: ifstream infile;
    infile.open(filename);
    if (!infile){
        std::cout<<"[ERROR] Couldn't open graph file\n";
        exit(1);
    } 

    int n,e;
    infile >> n >> e;

    this->__alloc__(n,e);
    for(int i=0;i<e;i++){
       infile >> this->edges[i];
    }

    for(int i=0;i<n;i++){
       infile >> this->offsets[i];
    }

    for(int i=0;i<e;i++){
       infile >> this->weights[i];
    }


}


template <class DataType, class FuncType >
void  Graph<DataType, FuncType> :: free()
{
    free(this->offsets);
    free(this->edges);
    free(this->weights);
}


#endif