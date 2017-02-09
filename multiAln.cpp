#include "AlnTree.h"


int main(int argc, char* argv[]){
	string seqfile = argv[1];
	string treefile = argv[2];
	float indel_rate = stof(argv[3]);
	string outfile = argv[4];
	
	SubMtrx m;
	AlnTree tree(m,indel_rate);
	cout << "Reading tree ... " << endl;
	tree.readNewick(treefile); 

	cout << "Reading sequences ... " << endl;
	tree.mapSeq2Tree(seqfile);

	cout << "Bottom up align ... " << endl;
	tree.btmupAln();

	cout << "Writing output ... " << endl;
	tree.printTopAln(outfile);
}
