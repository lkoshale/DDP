#include <iostream>
#include <cstdlib>
#include <string>

#include "graph.h"
#include "utils.h"

int main(){
    Graph<unsigned int > graph;
    std::string filename = "graph.txt";
    // graph.read_graph(filename);

    Builder<unsigned int> b;
    b.build_graph(build_fn);
    std::cout<<b.getN()<<" "<<b.getE()<<"\n";
    
    graph.build_graph(b);


    int* edge = graph.get_edges();
    int* off = graph.get_offsets();
    unsigned int* w = graph.get_weight();

    
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