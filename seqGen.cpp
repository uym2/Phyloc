#include "EvolTree.h"
#include <iostream>
using namespace std;

int main(){
	EvolTree tree;
	//tree.set_R("data/param.txt");
	//tree.print_R();
	tree.readNewick("data/avian-tree.tree");
	tree.writeNewick("out.tre");
	std::vector<double> pi = {0.292187, 0.211771, 0.229896, 0.266146};
	int N = 100000;
	int A = 0, C = 0, G = 0, T = 0;
	tree.seq_root(pi,N);
	for (int i = 0; i<tree.get_root()->get_seq().length(); i++){
		switch(tree.get_root()->get_seq()[i]){
			case 'A':
				A++;
				break;
			case 'C':
				C++;
				break;
			case 'G':
				G++;
				break;
			case 'T':
				T++;
				break;
		}
	}
	cout << (double)A/N << " " << (double)C/N << " " << (double)G/N << " " << (double)T/N;
	
	return 1;
}
