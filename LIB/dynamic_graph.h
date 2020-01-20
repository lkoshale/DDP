
#ifndef DYNAMIC_GRAPH_H
#define DYNAMIC_GRAPH_H

#include "graph.h"
#include "diff_graph.h"

template <class T = unsigned int>
class Dynamic_Graph
{
    private:
        Graph<T> *graph;
        Diff_Graph<T> *diff;
    
    public:
        Dynamic_Graph();

        Dynamic_Graph(Graph<T> *graph, Diff_Graph<T> * diff);

        Graph<T> get_graph(){return *(this->graph); }
        Diff_Graph<T> get_diff_graph(){ return *(this->diff); }

        void merge();
        void make_diff();
        void update_Del();

};


/*
* Implemetation
*/

template < class T > 
Dynamic_Graph<T> :: Dynamic_Graph()
{
    this->graph = new Graph();
    this->diff = new Diff_Graph(0);
}

template < class T >
Dynamic_Graph<T> :: Dynamic_Graph(Graph<T> *graph, Diff_Graph<T> * diff)
{
    this->graph = graph;
    this->diff = diff;
}

template < class T >
Graph<T> Dynamic_Graph<T> :: get_graph(){return *(this->graph); }

template < class T >
Diff_Graph<T> Dynamic_Graph<T> :: get_diff_graph(){ return *(this->diff); }



#endif