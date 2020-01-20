
#ifndef GRAPH_H
#define GRAPH_H


#include <unordered_map>


template <class T > 
class Graph 
{
    private:
        unsigned int N;         //nodes from 0 to N-1
        unsigned int E;
        int* offsets;
        int* edges;
        T* weights;

        std::unordered_map<long,int> vertex_map;     //if vertex doesn't follow 0-1 naming

    public: 
        Graph();

        void __alloc__(int n, int e);

        void build_graph(void (*builder)(int* offset,int* edges, T* weights) );

        int get_num_nodes();
        int get_num_edges();

        int* get_offsets();
        int* get_edges();
        T* get_weight();
        
};



/*
* Implementation Here as 
* templates cant have diff files for implementation
* (or move to another file .tpp nand include it here)
*/

template < class T > 
Graph<T>::Graph()
{
    N = 0;
    E = 0;
}

template < class T >
void Graph<T> :: __alloc__(int n, int e){
    N = n;
    E = e;
    this->offsets = (int*)malloc(sizeof(int)*N);
    this->edges = (int*)malloc(sizeof(int)*E);
    this->weights = (T*)malloc(sizeof(T)*E);
}

template < class T >
void Graph<T> :: build_graph(void (*builder)(int* offset,int* edges, T* weights) )
{
    builder(this->offsets,this->edges,this->weights);
}

template < class T >
int Graph<T> :: get_num_nodes(){ return N; }

template < class T >
int Graph<T> :: get_num_edges(){ return E; }

template < class T >
int* Graph<T> :: get_offsets(){ return this->offsets; }

template < class T >
int*  Graph<T> :: get_edges(){ return this->edges; }

template < class T >
T*  Graph<T> :: get_weight(){ return this->weights; }



#endif