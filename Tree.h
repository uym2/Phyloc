#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <fstream>
#include <boost/expm.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Aligner.h"

using namespace std;
using namespace boost::numeric::ublas;
class Node{
	std::vector<Node*> children;
	string label;
	double edge;
public:
	Node();
	Node(string label);
	~Node();
	
	bool set_label(string s);
	string get_label();
	
	double get_edge();
	bool set_edge(double e);

	Node* get_child(int n);	
	int nchild();

	bool writeNewick(ofstream &fout);
	bool add_child(Node *ch);
	bool is_leaf(){ return this->children.empty(); }
	
	virtual bool add_aln(string newAln) { return false;}
	
	virtual bool mapSeq(map<string,string> seq_map) { return false;}

	virtual float btmupAln(SubMtrx m, float idrate) { return false;}	

	virtual bool polyAln(SubMtrx m, float idrate) { return false;}	

	virtual string get_seq() { string s; return s;}
	
	virtual bool set_seq() { return false;}

	virtual bool transitAln(std::vector<string> &refAln, int shareIdx) { return false;}

	virtual bool add_gaps(int pos) { return false;}

	virtual bool printAln(ofstream &fout) { return false;}

	virtual std::vector<string>& getAln() {std::vector<string> foo; return foo; }
	
	virtual std::vector<string>& getLeafLabel() {std::vector<string> foo; return foo; }

	virtual bool set_seq(string seq) { return false;}
	
	virtual bool gen_seqs(const matrix<double> &R) { return false;}
	
	virtual char randit(const std::vector<double> &dtrb) { return 'x';};
};

class Tree{
	Node *root;
public:
	Tree();
	~Tree();
	bool readNewick(string treefile);
	bool writeNewick(string treefile);
	Node*& get_root() { return this->root; }
	virtual Node* create_node() = 0;
};
