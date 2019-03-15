#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <queue>


using namespace std;


typedef struct node {
	int val;
	struct node* edge[10000];
	int size;
	bool visited;
	int level;

}Node;


Node* create_node(int val){
	Node* temp = (Node*)malloc(sizeof(Node));
	if(temp!=NULL){
		temp->val = val;
		temp->size = 0;
		temp->visited = false;
		temp->level = -1;
	}

	return temp;
}

/* depth first traversal */
void dfs(Node* root){
	if(root==NULL)
		return;

	printf("%d\n",root->val);
	root->visited = true;
	for(int i=0;i<root->size;i++){
		if(root->edge[i]->visited==false)
			dfs(root->edge[i]);
	}
}


void create_edge(Node* a,Node* b){
	a->edge[a->size] = b;
	a->size++;
	// b->edge[b->size] = a;
	// b->size++;
}

void bfs(Node* root){
	if(root==NULL)
		return;

	root->level = 0;
	queue<Node*> qu;
	qu.push(root);
	root->visited=true;
	while(!qu.empty()){
		Node* t = qu.front();
		qu.pop();
		// printf("%d\n",t->val);
		for(int i=0;i<t->size;i++){
			if(t->edge[i]->visited==false){
				qu.push(t->edge[i]);
				t->edge[i]->visited = true;
				t->edge[i]->level = t->level +1;
			}
		}
	}
}

int main(int argc, char const *argv[])
{
	int N;
	cin>>N;
	int startNode;
	cin>>startNode;

	vector<Node*> vec;
	for(int i=0;i<N;i++){
		vec.push_back(create_node(i));
	}

	int a,b;
	while(scanf("%d %d\n",&a,&b)!=EOF){
		if(a>=N ||b>= N ){
			printf("invalid input\n");
			printf("%d %d\n",a,b);
			exit(0);
		}
		create_edge(vec[a],vec[b]);
	}

	// dfs(vec[0]);
	// for(int i=0;i<vec.size();i++){
	// 	vec[i]->visited = false;
	// }
	// printf("\n\n");
	bfs(vec[startNode]);
	int maxL = 0;
	for(int i=0;i<N;i++){
		if( vec[i]->level > maxL)
			maxL = vec[i]->level;

		// cout<<i<<" "<<vec[i]->level<<"\n";
	} 
	cout<<maxL<<"\n";
	return 0;
}