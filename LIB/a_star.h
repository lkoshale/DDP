#ifndef A_STAR_H
#define A_STAR_H

#include "dynamic_graph.h"

template <class T, class U >
class A_star
{
    private:
        Dynamic_graph<T> *graph;
        Dynamic_graph<T> rev_graph;
        U* Hx;
        U* Cx;
        unsigned int* PQ;
        unsigned int num_pq;
        unsigned int start;
        unsigned int end;
    
    public:
        A_star(Dynamic_graph<T> *graph, unsigned int start,unsigned int end, unsigned int K ){
            this->graph = graph;
            this->num_pq = K;
            this->start = start;
            this->end = end;
        }

        std::vector<unsigned int> get_path_gpu();
    
};

#endif