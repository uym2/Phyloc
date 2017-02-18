#include "EvolTree.h"
#include <iostream>
using namespace std;

int main(){
	EvolTree tree;
	tree.set_R("data/param.txt");
	tree.print_R();
	//tree.readNewick("data/avian-tree.tree");
	//tree.writeNewick("out.tre");

	return 1;
}
