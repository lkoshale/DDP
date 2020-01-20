#ifdef __NVCC__

#ifndef DYNAMIC_GRAPH_CUH
#define DYNAMIC_GRAPH_CUH

#include "diff_graph.cuh"

template < class T >
class GPU_Dynamic_Graph
{
    private:
        GPU_Graph<T> *graph;
        GPU_Diff_Graph<T> *diff;

    public:
        GPU_Dynamic_Graph();

        GPU_Dynamic_Graph(GPU_Graph<T> * graph, GPU_Diff_Graph<T> * diff);

        ~GPU_Diff_Graph();

};



template < class T >
GPU_Dynamic_Graph<T> :: GPU_Dynamic_Graph()
{
    this->graph = new GPU_Graph();
    this->diff = new GPU_Diff_Graph(0);
}


#endif

#endif