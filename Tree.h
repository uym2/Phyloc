#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <fstream>
#include "Aligner.h"

using namespace std;

class Node{
	typedef Node myNode;
	vector<myNode*> children;
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

	myNode*& get_child(int n) { return this->children[n]; }	
	int nchild();

	bool writeNewick(ofstream &fout);
	bool add_child(myNode *ch);
	bool is_leaf(){ return this->children.empty(); }
	
	virtual bool add_aln(string newAln) { return false; }
	
	virtual bool mapSeq(map<string,string> seq_map)  { return false; }

	virtual float btmupAln(SubMtrx m, float idrate)  { return false; }	

	virtual bool polyAln(SubMtrx m, float idrate)  { return false; }

	virtual string get_seq() { return "";}
	
	virtual bool set_seq() { return false; }

	virtual bool transitAln(vector<string> &refAln, int shareIdx) { return false; }

	virtual bool add_gaps(int pos) {return false; }

	virtual bool printAln(ofstream &fout) { return false; }

	virtual vector<string>& getAln() { vector<string> foo; return foo; }
	
	virtual vector<string>& getLeafLabel() {vector<string> foo; return foo; }
};

class Tree{
	typedef Node myNode;
	myNode *root;
public:
	Tree();
	~Tree();
	bool readNewick(string treefile);
	bool writeNewick(string treefile);
	myNode*& get_root() { return this->root; }
	virtual myNode* create_node() = 0;
};
