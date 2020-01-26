#include <iostream>
#include <cstdlib>
#include <string>

#include "graph.h"
#include "diff_graph.h"

#ifdef __NVCC__
    #include <cuda.h>
    #include "graph.cuh"
#endif


int main(){
    Graph<unsigned int > graph;
    std::string filename = "graph.txt";
    graph.read_graph(filename);
    
   
    int* edge = graph.get_edges();
    int* off = graph.get_offsets();
    unsigned int* w = graph.get_weights();

    #ifdef __NVCC__
        GPU_Graph<unsigned int> g(&graph);
        printf("allocated\n");
    #endif
 
    
    for(int i=0;i<graph.get_num_edges();i++){
        std::cout<<edge[i]<<" ";
    }
    std::cout<<"\n";

    for(int i=0;i<graph.get_num_nodes();i++){
        std::cout<<off[i]<<" ";
    }
    std::cout<<"\n";

    for(int i=0;i<graph.get_num_edges();i++){
        std::cout<<w[i]<<" ";
    }
    std::cout<<"\n";
    

    return 0;
}