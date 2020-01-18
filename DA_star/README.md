
## Issue 
 - contains the logical problems in implementation.

## outline.txt
 - outline of how DA_insert.cu works.


## genDynamicGraph.cpp 
 - takes graph in form u v w
 - removes some of the edges as for insert
 - mark some of the edges to be deleted
 - outputs the update.txt and  graph.txt ( graph with edges to be added in update)


## genOpCSR.cpp 
 - takes the graph.txt and genrates the garph where edges are reversed 
 - outputs graph_op.txt

### genCSRwholegraph.cpp
 - genrates graph with all deletions from Updates.txt and all insertions added.

## bfs_hx.cu 
 - takes graph.txt and do bfs from input startNode
 - outputs Hx.txt ( used for doing bfs on back edge garph to generate Hx)

## DA_v0.cu
 - Does A* on graph.txt and genrated diff graph from Updates.txt
 - propgate the updates and does A* again 


#DA_dyn_final.cu
	- single propagation

#DA_dyn_cpu.cu
	- separate propgation
	- no cycle
	- on cpu
	- working
	
#DA_dyn_ins_del.cu
	- CX defined as global var
	- canot copy so cant do cycle check on cpu

#DA_dyn_ins_del_nocycle.cu
	- working
	- moved everything onto gpu


	
#DA_dyn_ins_del.cu
