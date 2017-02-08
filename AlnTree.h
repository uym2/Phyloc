#include "Tree.h"

class AlnNode: public Node{
	vector<string> Aln;
public:
	AlnNode(){
		this->Aln = NULL;
	}
	~AlnNode(){
	}
};

class AlnTree: public Tree{
public:
	Node* create_node(){
		Node* p = new AlnNode;
		return p;
	}
};
