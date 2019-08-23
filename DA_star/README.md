
# genDynamicGraph.cpp 
 - takes graph in form u v w
 - removes some of the edges as for insert
 - outputs the update.txt and  graph.txt ( graph with edges to be added in update)


# genOpCSR.cpp 
 - takes the graph.txt and genrates the garph where edges are reversed 
 - outputs graph_op.txt


# bfs_hx.cu 
 - takes graph.txt and do bfs from input startNode
 - outputs Hx.txt ( used for doing bfs on back edge garph to generate Hx)

# DA_v0.cu
 - Does A* on graph.txt and genrated diff graph from Updates.txt
 - propgate the updates and does A* again 
 