#ifdef __NVCC__

#ifndef GRAPH_CUH
#define GRAPH_CUH

#include "utils.cuh"
#include "graph.h"


template < class T = unsigned int>
class GPU_Graph : public Graph
{
    private:
        int* d_offsets;
        int* d_edges;
        T* d_weights;
    
    public:

        GPU_Graph();
         
        void __gpu_alloc();

        void __gpu_free();

        int* get_offsets(){ return this->d_offsets; }
        int* get_edges(){ return this->d_edges; }
        T* get_weight(){ return this->d_weights; }

};


/*
* Implemenattion
* Here
*/


template< class T >
GPU_Graph<T>:: GPU_Graph()
{
    Base::Graph();
}

template< class T >
void GPU_Graph<T>:: __gpu_alloc()
{
    gpuErrchk ( cudaMalloc(&d_offsets,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&d_edges,sizeof(int)*E) );
    gpuErrchk ( cudaMalloc(&d_weights,sizeof(T)*E) );

    gpuErrchk ( cudaMemcpy(d_offsets,offsets,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_edges,edges,sizeof(int)*E,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_weights,weights,sizeof(T)*E,cudaMemcpyHostToDevice) );

}

template< class T >
int* GPU_Graph<T>:: get_offsets(){ return this->d_offsets; }

template< class T >
int* GPU_Graph<T>:: get_edges(){ return this->d_edges; }

template< class T >
T* GPU_Graph<T>:: get_weight(){ return this->d_weights; }


template< class T >
void  GPU_Graph<T>:: __gpu_free(){
    gpuErrchk ( cudaFree(d_offsets) );
    gpuErrchk ( cudaFree(d_edges) );
    gpuErrchk ( cudaFree(d_weights) );
}



#endif

#endif
