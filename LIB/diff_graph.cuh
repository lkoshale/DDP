
#ifdef __NVCC__

#ifndef DIFF_GRAPH_CUH
#define DIFF_GRAPH_CUH

#include "graph.h"

template < class T = unsigned int>
class GPU_Diff_Graph : GPU_Graph
{
    private:
        int id;
    
    public:
        GPU_Diff_Graph(int id);
};


/*************************
*       Implementation
**************************/

template< class T >
GPU_Diff_Graph<T> :: GPU_Diff_Graph(int id){
    this->id = id;
    GPU_Graph();
}


#endif

#endif